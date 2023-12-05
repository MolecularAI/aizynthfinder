import pytest

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.policy import (
    MultiExpansionStrategy,
    TemplateBasedExpansionStrategy,
)
from aizynthfinder.utils.exceptions import PolicyException


def test_multi_expansion_strategy_incorrect_keys(
    default_config, setup_template_expansion_policy
):
    expansion_policy = default_config.expansion_policy
    strategy1, _ = setup_template_expansion_policy("policy1")
    expansion_policy.load(strategy1)

    with pytest.raises(ValueError, match="expansion strategy keys must exist"):
        multi_expansion_strategy = MultiExpansionStrategy(
            "multi_expansion_strategy",
            default_config,
            expansion_strategies=["policy1", "policy2"],
        )
        mols = [TreeMolecule(smiles="CCO", parent=None)]
        multi_expansion_strategy.get_actions(mols)


def test_multi_expansion_strategy(default_config, setup_template_expansion_policy):
    expansion_policy = default_config.expansion_policy
    strategy1, _ = setup_template_expansion_policy("policy1")
    expansion_policy.load(strategy1)
    strategy2, _ = setup_template_expansion_policy("policy2")
    expansion_policy.load(strategy2)
    strategy3, _ = setup_template_expansion_policy("policy3")
    expansion_policy.load(strategy3)

    multi_expansion_strategy = MultiExpansionStrategy(
        "multi_expansion_strategy",
        default_config,
        expansion_strategies=["policy1", "policy2"],
    )
    multi_expansion_strategy.additive_expansion = True

    mols = [TreeMolecule(smiles="CCO", parent=None)]
    _, priors = multi_expansion_strategy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.7, 0.2]


def test_multi_expansion_strategy_wo_additive_expansion(
    default_config, setup_template_expansion_policy
):
    expansion_policy = default_config.expansion_policy
    strategy1, _ = setup_template_expansion_policy("policy1")
    expansion_policy.load(strategy1)
    strategy2, _ = setup_template_expansion_policy("policy2")
    expansion_policy.load(strategy2)

    multi_expansion_strategy = MultiExpansionStrategy(
        "multi_expansion_strategy",
        default_config,
        expansion_strategies=["policy1", "policy2"],
    )

    mols = [TreeMolecule(smiles="CCO", parent=None)]
    _, priors = multi_expansion_strategy.get_actions(mols)

    assert priors == [0.7, 0.2]


def test_create_templated_expansion_strategy_wo_kwargs():
    with pytest.raises(
        PolicyException, match=" class needs to be initiated with keyword arguments"
    ):
        _ = TemplateBasedExpansionStrategy("dummy", None)


def test_load_templated_expansion_strategy(
    default_config, setup_template_expansion_policy
):
    strategy, mocked_onnx_model = setup_template_expansion_policy()
    mocked_onnx_model.assert_called_once()
    assert len(strategy.templates) == 3


def test_load_invalid_templated_expansion_strategy(
    default_config, create_dummy_templates, mock_onnx_model
):
    templates_filename = create_dummy_templates(4)
    with pytest.raises(PolicyException):
        TemplateBasedExpansionStrategy(
            "policy1",
            default_config,
            model="dummy.onnx",
            template=templates_filename,
        )


def test_load_templated_expansion_strategy_from_csv(
    default_config, mock_onnx_model, tmpdir
):
    templates_filename = str(tmpdir / "temp.csv")

    with open(templates_filename, "w") as fileobj:
        fileobj.write("template_index\ttemplate\tmetadata\n")
        fileobj.write("0\tAAA\tmetadata1\n")
        fileobj.write("1\tBBB\tmetadata2\n")
        fileobj.write("2\tCCC\tmetadata3\n")

    strategy = TemplateBasedExpansionStrategy(
        "default", default_config, model="dummy.onnx", template=templates_filename
    )

    assert len(strategy.templates) == 3
    assert list(strategy.templates.columns) == ["template", "metadata"]
