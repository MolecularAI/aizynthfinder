import pytest

from aizynthfinder.chem import Molecule
from aizynthfinder.mcts.policy import PolicyException


def test_load_policy(policy, mocker, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    mocked_load_model = mocker.patch("aizynthfinder.utils.keras_utils.load_model")

    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    mocked_load_model.assert_called_once_with("dummy.hdf5", custom_objects=mocker.ANY)
    assert len(policy["policy1"]["templates"]) == 3


def test_select_policy(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    mocked_model = mock_policy_model
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    policy.select_policy("policy1")

    assert policy._policy_model.model is mocked_model.return_value
    assert len(policy._templates) == 3


def test_select_policy_invalid_key(policy, mock_policy_model, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")

    with pytest.raises(PolicyException):
        policy.select_policy("policy2")


def test_get_actions(policy, mock_policy_model, mock_stock, mocker, shared_datadir):
    templates_filename = str(shared_datadir / "templates.hdf5")
    policy.load_policy("dummy.hdf5", templates_filename, "policy1")
    policy.select_policy("policy1")
    mols = [Molecule(smiles="CCO")]
    mocker.patch("aizynthfinder.mcts.stock.Stock.__contains__", return_value=False)

    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.2]

    policy._config.cutoff_cumulative = 1.0
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7, 0.2, 0.1]

    policy._config.cutoff_number = 1
    actions, priors = policy.get_actions(mols)

    assert priors == [0.7]
