import pytest
import numpy as np

from aizynthfinder.chem import Molecule, TreeMolecule
from aizynthfinder.context.policy import PolicyException


def test_load_policy(expansion_policy, mocker, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")

    expansion_policy.load("dummy.hdf5", templates_filename, "policy1")

    mock_policy_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)
    assert len(expansion_policy["policy1"]["templates"]) == 3

    # Now load something with a different output dimension
    mock_policy_model.return_value.output = np.zeros((2, 2))
    with pytest.raises(PolicyException):
        expansion_policy.load("dummy.hdf5", templates_filename, "policy1")


def test_get_actions(
    expansion_policy, mock_policy_model, mock_stock, mocker, shared_datadir
):
    templates_filename = str(shared_datadir / "templates.hdf5")
    expansion_policy.load("dummy.hdf5", templates_filename, "policy1")
    expansion_policy.select("policy1")
    mols = [Molecule(smiles="CCO")]
    mocker.patch("aizynthfinder.context.stock.Stock.__contains__", return_value=False)

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


def test_get_actions_two_policies(
    expansion_policy, mock_policy_model, mocker, shared_datadir
):
    templates_filename = str(shared_datadir / "templates.hdf5")
    expansion_policy.load("dummy.hdf5", templates_filename, "policy1")
    expansion_policy.load("dummy.hdf5", templates_filename, "policy2")
    expansion_policy.select(["policy1", "policy2"])
    mols = [Molecule(smiles="CCO")]
    mocker.patch("aizynthfinder.context.stock.Stock.__contains__", return_value=False)

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


def test_load_filter_policy(filter_policy, mocker):
    mocked_load_model = mocker.patch("aizynthfinder.utils.models.load_keras_model")

    filter_policy.load("dummy.hdf5", "policy1")

    mocked_load_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)


def test_feasible(filter_policy, mock_policy_model, mocker, simple_actions):
    filter_policy.load("dummy.hdf5", "policy1")
    filter_policy.select("policy1")
    mol = TreeMolecule(parent=None, smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    reactions, _ = simple_actions(mol)

    filter_policy._config.filter_cutoff = 0.9
    assert not filter_policy.is_feasible(reactions[0])

    filter_policy._config.filter_cutoff = 0.15
    assert filter_policy.is_feasible(reactions[0])
