import pytest
import numpy as np

from aizynthfinder.chem import TreeMolecule, SmilesBasedRetroReaction
from aizynthfinder.context.policy import (
    TemplateBasedExpansionStrategy,
    QuickKerasFilter,
)
from aizynthfinder.utils.exceptions import RejectionException, PolicyException


def test_create_templated_expansion_strategy_wo_kwargs():
    with pytest.raises(
        PolicyException, match=" class needs to be initiated with keyword arguments"
    ):
        _ = TemplateBasedExpansionStrategy("dummy", None)


def test_load_templated_expansion_policy(
    default_config, setup_template_expansion_policy, mocker
):
    strategy, mocked_keras_model = setup_template_expansion_policy()
    mocked_keras_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)
    assert len(strategy.templates) == 3


def test_load_invalid_templated_expansion_policy(
    default_config, create_dummy_templates, mock_keras_model
):
    templates_filename = create_dummy_templates(3)
    mock_keras_model.return_value.output = np.zeros((2, 2))
    with pytest.raises(PolicyException):
        TemplateBasedExpansionStrategy(
            "policy1",
            default_config,
            source="dummy.hdf5",
            templatefile=templates_filename,
        )


def test_load_expansion_policy(default_config, setup_template_expansion_policy):
    strategy, _ = setup_template_expansion_policy()
    expansion_policy = default_config.expansion_policy
    expansion_policy.load(strategy)

    with pytest.raises(PolicyException):
        expansion_policy.load(5)


def test_load_expansion_policy_from_config_files(
    default_config, mock_keras_model, create_dummy_templates
):
    template_filename = create_dummy_templates(3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load_from_config(
        **{
            "files": {
                "policy1": ["dummy1", template_filename],
                "policy2": ["dummy1", template_filename],
            }
        }
    )
    assert "policy1" in expansion_policy.items
    assert len(expansion_policy["policy1"].templates) == 3
    assert "policy2" in expansion_policy.items
    assert len(expansion_policy["policy2"].templates) == 3


def test_load_expansion_policy_from_config_custom(
    default_config, mock_keras_model, create_dummy_templates
):
    template_filename = create_dummy_templates(3)
    expansion_policy = default_config.expansion_policy
    expansion_policy.load_from_config(
        **{
            "TemplateBasedExpansionStrategy": {
                "policy1": {"source": "dummy1", "templatefile": template_filename}
            },
            "aizynthfinder.context.policy.TemplateBasedExpansionStrategy": {
                "policy2": {"source": "dummy1", "templatefile": template_filename}
            },
        }
    )
    assert "policy1" in expansion_policy.items
    assert len(expansion_policy["policy1"].templates) == 3
    assert "policy2" in expansion_policy.items
    assert len(expansion_policy["policy2"].templates) == 3


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

    expansion_policy._config.cutoff_cumulative = 1.0
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1]

    expansion_policy._config.cutoff_number = 1
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7]


def test_get_actions_two_policies(default_config, setup_template_expansion_policy):
    expansion_policy = default_config.expansion_policy
    strategy1, _ = setup_template_expansion_policy("policy1")
    expansion_policy.load(strategy1)
    strategy2, _ = setup_template_expansion_policy("policy2")
    expansion_policy.load(strategy2)
    default_config.additive_expansion = True

    expansion_policy.select(["policy1", "policy2"])
    mols = [TreeMolecule(smiles="CCO", parent=None)]

    actions, priors = expansion_policy.get_actions(mols)

    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1"] * 2 + ["policy2"] * 2
    assert priors == [0.7, 0.2, 0.7, 0.2]

    expansion_policy._config.cutoff_cumulative = 1.0
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1, 0.7, 0.2, 0.1]

    expansion_policy._config.cutoff_number = 1
    actions, priors = expansion_policy.get_actions(mols)

    assert priors == [0.7, 0.7]

    default_config.additive_expansion = False
    default_config.cutoff_number = 2

    actions, priors = expansion_policy.get_actions(mols)

    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1", "policy1"]
    assert priors == [0.7, 0.2]


def test_create_quick_filter_strategy_wo_kwargs():
    with pytest.raises(
        PolicyException, match=" class needs to be initiated with keyword arguments"
    ):
        _ = QuickKerasFilter("dummy", None)


def test_load_filter_policy(default_config, mock_keras_model, mocker):
    strategy = QuickKerasFilter("policy1", default_config, source="dummy.hdf5")
    default_config.filter_policy.load(strategy)

    mock_keras_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)

    with pytest.raises(PolicyException):
        default_config.filter_policy.load(5.0)


def test_load_filter_policy_from_config_files(default_config, mock_keras_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{
            "files": {
                "policy1": "dummy1",
                "policy2": "dummy1",
            }
        }
    )
    assert "policy1" in filter_policy.items
    assert "policy2" in filter_policy.items


def test_load_filter_policy_from_config_custom(default_config, mock_keras_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(
        **{
            "QuickKerasFilter": {"policy1": {"source": "dummy1"}},
            "aizynthfinder.context.policy.QuickKerasFilter": {
                "policy2": {"source": "dummy1"}
            },
        }
    )
    assert "policy1" in filter_policy.items
    assert "policy2" in filter_policy.items


def test_filter_rejection(default_config, mock_keras_model):
    filter_policy = default_config.filter_policy
    filter_policy.load_from_config(**{"files": {"policy1": "dummy1"}})
    mol = TreeMolecule(
        parent=None, smiles="CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    )
    reaction = SmilesBasedRetroReaction(
        mol, reactants_str="CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O"
    )

    with pytest.raises(PolicyException, match="selected"):
        filter_policy(reaction)

    filter_policy.select("policy1")
    filter_policy._config.filter_cutoff = 0.9
    with pytest.raises(RejectionException):
        filter_policy(reaction)

    filter_policy._config.filter_cutoff = 0.15
    filter_policy(reaction)
