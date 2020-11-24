import pytest
import numpy as np

from aizynthfinder.chem import Molecule
from aizynthfinder.mcts.policy import PolicyException


def test_load_policy(policy, mocker, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")

    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    mock_policy_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)
    assert len(policy["policy1"]["templates"]) == 3

    # Now load something with a different output dimension
    mock_policy_model.return_value.output = np.zeros((2, 2))
    with pytest.raises(PolicyException):
        policy.load_policy("dummy.hdf5", templates_filename, "policy1")


def test_select_policy(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    policy.select_policy("policy1")

    assert len(policy.selected_policies) == 1
    assert policy.selected_policy == "policy1"


def test_select_policies_only_one(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    policy.select_policies("policy1")

    assert len(policy.selected_policies) == 1


def test_select_two_policies(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")
    policy.load_policy("dummy.hdf5", templates_filename, "policy2")

    policy.select_policies(["policy1", "policy2"])

    assert len(policy.selected_policies) == 2


def test_select_policy_invalid_key(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    with pytest.raises(PolicyException):
        policy.select_policy("policy2")

    with pytest.raises(PolicyException):
        policy.select_policies("policy2")

    with pytest.raises(PolicyException):
        policy.select_policies(["policy2", "policy3"])


def test_get_actions(policy, mock_policy_model, mocker, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")
    policy.select_policy("policy1")
    mols = [Molecule(smiles="CCO")]
    mocker.patch("aizynthfinder.mcts.stock.Stock.__contains__", return_value=False)

    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.2]
    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1", "policy1"]

    policy._config.cutoff_cumulative = 1.0
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1]

    policy._config.cutoff_number = 1
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7]


def test_get_actions_two_policies(policy, mock_policy_model, mocker, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")
    policy.load_policy("dummy.hdf5", templates_filename, "policy2")
    policy.select_policies(["policy1", "policy2"])
    mols = [Molecule(smiles="CCO")]
    mocker.patch("aizynthfinder.mcts.stock.Stock.__contains__", return_value=False)

    actions, priors = policy.get_actions(mols)

    policy_names = [action.metadata["policy_name"] for action in actions]
    assert policy_names == ["policy1"] * 2 + ["policy2"] * 2
    assert priors == [0.7, 0.2, 0.7, 0.2]

    policy._config.cutoff_cumulative = 1.0
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1, 0.7, 0.2, 0.1]

    policy._config.cutoff_number = 1
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.7]
