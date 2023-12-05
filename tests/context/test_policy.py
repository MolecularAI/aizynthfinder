import numpy as np
import pytest

from aizynthfinder.chem import (
    SmilesBasedRetroReaction,
    TemplatedRetroReaction,
    TreeMolecule,
)
from aizynthfinder.context.policy import (
    QuickKerasFilter,
    ReactantsCountFilter,
)
from aizynthfinder.utils.exceptions import PolicyException, RejectionException


def test_load_expansion_policy(default_config, setup_template_expansion_policy):
    strategy, _ = setup_template_expansion_policy()
    expansion_policy = default_config.expansion_policy
    expansion_policy.load(strategy)

    with pytest.raises(PolicyException):
        expansion_policy.load(5)


def test_load_expansion_policy_from_config_files(
    default_config, mock_onnx_model, create_dummy_templates
):
    template_filename = create_dummy_templates(3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load_from_config(
        **{
            "policy1": {
                "type": "template-based",
                "model": "dummy1.onnx",
                "template": template_filename,
            },
            "policy2": {
                "type": "template-based",
                "model": "dummy1.onnx",
                "template": template_filename,
            },
        }
    )
    assert "policy1" in expansion_policy.items
    assert len(expansion_policy["policy1"].templates) == 3
    assert "policy2" in expansion_policy.items
    assert len(expansion_policy["policy2"].templates) == 3


def test_load_expansion_policy_from_config_custom(
    default_config, mock_onnx_model, create_dummy_templates
):
    template_filename = create_dummy_templates(3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load_from_config(
        **{
            "policy1": {
                "type": "TemplateBasedExpansionStrategy",
                "model": "dummy1.onnx",
                "template": template_filename,
                "cutoff_number": 75,
            },
            "policy2": {
                "type": "aizynthfinder.context.policy.TemplateBasedExpansionStrategy",
                "model": "dummy1.onnx",
                "template": template_filename,
                "cutoff_number": 25,
            },
            "policy3": {
                "type": "aizynthfinder.context.policy.TemplateBasedExpansionStrategy",
                "model": "dummy1.onnx",
                "template": template_filename,
            },
        }
    )
    assert "policy1" in expansion_policy.items
    assert len(expansion_policy["policy1"].templates) == 3
    assert "policy2" in expansion_policy.items
    assert len(expansion_policy["policy2"].templates) == 3
    assert "policy3" in expansion_policy.items
    assert len(expansion_policy["policy3"].templates) == 3
    assert expansion_policy["policy1"].cutoff_number == 75
    assert expansion_policy["policy2"].cutoff_number == 25
    assert expansion_policy["policy3"].cutoff_number == 50


def test_get_actions(default_config, setup_template_expansion_policy):
    strategy, _ = setup_template_expansion_policy()
    expansion_policy = default_config.expansion_policy
    expansion_policy.load(strategy)

    mols = [TreeMolecule(smiles="CCO", parent=None)]

    with pytest.raises(PolicyException, match="selected"):
        expansion_policy.get_actions(mols)

    expansion_policy.select("policy1")
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.2]
    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1", "policy1"]

    default_config.expansion_policy["policy1"].cutoff_cumulative = 1.0
    expansion_policy.reset_cache()
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1]

    default_config.expansion_policy["policy1"].cutoff_number = 1
    expansion_policy.reset_cache()
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7]

    default_config.expansion_policy["policy1"].cutoff_number = 2
    default_config.expansion_policy["policy1"].rescale_prior = True
    expansion_policy.reset_cache()
    actions, priors = expansion_policy.get_actions(mols)
    assert [round(prior, 1) for prior in priors] == [0.6, 0.4]


def test_get_actions_two_policies(default_config, setup_template_expansion_policy):
    expansion_policy = default_config.expansion_policy
    strategy1, _ = setup_template_expansion_policy("policy1")
    expansion_policy.load(strategy1)
    strategy2, _ = setup_template_expansion_policy("policy2")
    expansion_policy.load(strategy2)

    expansion_policy.select(["policy1", "policy2"])
    mols = [TreeMolecule(smiles="CCO", parent=None)]

    actions, priors = expansion_policy.get_actions(mols)

    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1"] * 2 + ["policy2"] * 2
    assert priors == [0.7, 0.2, 0.7, 0.2]

    expansion_policy["policy1"].cutoff_cumulative = 1.0
    expansion_policy["policy2"].cutoff_cumulative = 1.0
    expansion_policy.reset_cache()
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1, 0.7, 0.2, 0.1]

    expansion_policy["policy1"].cutoff_number = 1
    expansion_policy["policy2"].cutoff_number = 1
    expansion_policy.reset_cache()
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.7]


def test_get_actions_using_rdkit(
    default_config, setup_template_expansion_policy, mocker
):
    smarts = (
        "([#8:4]-[N;H0;D3;+0:5](-[C;D1;H3:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])"
        ">>(Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([#8:4]-[NH;D2;+0:5]-[C;D1;H3:6])"
    )
    strategy, _ = setup_template_expansion_policy(templates=[smarts] * 3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load(strategy)
    mols = [TreeMolecule(smiles="CCO", parent=None)]
    expansion_policy.select("policy1")
    mocker.patch(
        "aizynthfinder.chem.reaction.AllChem.ReactionFromSmarts",
        side_effect=RuntimeError("Intential error"),
    )

    actions, _ = expansion_policy.get_actions(mols)

    assert actions[0] is not None

    # Now switch to RDKit
    default_config.expansion_policy["policy1"].use_rdchiral = False
    actions, _ = expansion_policy.get_actions(mols)

    # This is a kind of convulted way to check that the RDKit application route will be called
    # the actual testing of the RDKit routine is done elsewhere
    with pytest.raises(RuntimeError, match="Intential error"):
        _ = actions[0].reactants


def test_template_based_expansion_caching(
    default_config, mock_onnx_model, create_dummy_templates
):
    template_filename = create_dummy_templates(3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load_from_config(
        **{"policy1": {"model": "dummy1.onnx", "template": template_filename}},
    )
    policy = expansion_policy["policy1"]
    mols = [TreeMolecule(smiles="CCO", parent=None)]
    caching_mols = [TreeMolecule(smiles="CCCCO", parent=None)]

    actions1, priors1 = policy(mols, caching_mols)
    actions2, priors2 = policy(mols)

    mock_onnx_model.assert_called_once()
    assert priors1 == priors2
    assert actions1[0].smarts == actions2[0].smarts


def test_create_quick_filter_strategy_wo_kwargs():
    with pytest.raises(
        PolicyException, match=" class needs to be initiated with keyword arguments"
    ):
        _ = QuickKerasFilter("dummy", None)


def test_load_filter_policy(default_config, mock_onnx_model):
    strategy = QuickKerasFilter("policy1", default_config, model="dummy.onnx")
    default_config.filter_policy.load(strategy)

    mock_onnx_model.assert_called_once()

    with pytest.raises(PolicyException):
        default_config.filter_policy.load(5.0)


def test_load_filter_policy_from_config_files(default_config, mock_onnx_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{
            "policy1": {"type": "quick-filter", "model": "dummy1.onnx"},
            "policy2": {"type": "quick-filter", "model": "dummy1.onnx"},
            "policy3": "dummy1.onnx",
        }
    )
    assert "policy1" in filter_policy.items
    assert "policy2" in filter_policy.items
    assert "policy3" in filter_policy.items


def test_load_filter_policy_from_config_custom(default_config, mock_onnx_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{
            "policy1": {"type": "QuickKerasFilter", "model": "dummy1.onnx"},
            "policy2": {"type": "feasibility", "model": "dummy1.onnx"},
        }
    )
    assert "policy1" in filter_policy.items
    assert "policy2" in filter_policy.items


def test_filter_rejection(default_config, mock_onnx_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{"policy1": {"type": "quick-filter", "model": "dummy1.onnx"}}
    )
    mol = TreeMolecule(
        parent=None, smiles="CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    )
    reaction = SmilesBasedRetroReaction(
        mol, reactants_str="CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O"
    )

    with pytest.raises(PolicyException, match="selected"):
        filter_policy(reaction)

    filter_policy.select("policy1")
    filter_policy._config.filter_policy["policy1"].filter_cutoff = 0.9
    with pytest.raises(RejectionException):
        filter_policy(reaction)

    filter_policy._config.filter_policy["policy1"].filter_cutoff = 0.15
    filter_policy(reaction)


def test_skip_filter_rejection(default_config, mock_onnx_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{
            "policy1": {
                "type": "quick-filter",
                "model": "dummy1.onnx",
                "exclude_from_policy": ["dummy.onnx"],
            },
        }
    )
    mol = TreeMolecule(
        parent=None, smiles="CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    )
    filter_policy.select("policy1")
    filter_policy._config.filter_policy["policy1"].filter_cutoff = 0.9

    reaction = SmilesBasedRetroReaction(
        mol, reactants_str="CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O"
    )

    with pytest.raises(RejectionException):
        filter_policy(reaction)

    reaction = SmilesBasedRetroReaction(
        mol,
        reactants_str="CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O",
        metadata={"policy_name": "dummy.onnx"},
    )
    assert filter_policy(reaction) is None


def test_reactants_count_rejection(default_config):
    smarts = (
        "([C:3]-[N;H0;D2;+0:2]=[C;H0;D3;+0:1](-[c:4]1:[c:5]:[c:6]:[c:7]:[c:8]:[c:9]:1)-[c;H0;D3;+0:11](:[c:10]):[c:12])>>"
        "(O=[C;H0;D3;+0:1](-[NH;D2;+0:2]-[C:3])-[c:4]1:[c:5]:[c:6]:[c:7]:[c:8]:[c:9]:1.[c:10]:[cH;D2;+0:11]:[c:12])"
    )
    mol = TreeMolecule(parent=None, smiles="c1c2c(ccc1)CCN=C2c3ccccc3")
    rxn1 = TemplatedRetroReaction(mol=mol, smarts=smarts)
    filter = ReactantsCountFilter("dummy", default_config)

    assert len(rxn1.reactants) == 2

    rxn2 = rxn1.copy(index=1)
    if len(rxn1.reactants[0]) == 1:
        rxn1, rxn2 = rxn2, rxn1

    assert filter(rxn2) is None

    with pytest.raises(RejectionException):
        filter(rxn1)
