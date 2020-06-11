import sys

import pytest
import numpy as np
import yaml

from aizynthfinder.mcts.config import Configuration
from aizynthfinder.mcts.policy import Policy
from aizynthfinder.mcts.stock import Stock
from aizynthfinder.mcts.node import Node
from aizynthfinder.chem import Reaction
from aizynthfinder.mcts.mcts import SearchTree


def pytest_addoption(parser):
    parser.addoption(
        "--finder_config", help="the configuration file for the aizynthfinder",
    )
    parser.addoption(
        "--stocks", nargs="+", help="the stocks to use in the aizynthfinder",
    )
    parser.addoption("--policy", help="the policy to use in the aizynthfinder")
    parser.addoption(
        "--run_integration",
        action="store_true",
        default=False,
        help="run integration tests",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run_integration"):
        return

    skip_integration = pytest.mark.skip(reason="need --run_integration option to run")
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(skip_integration)


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "integration: this one is for integration tests."
    )


@pytest.fixture
def add_cli_arguments():
    saved_argv = list(sys.argv)

    def wrapper(args):
        sys.argv = [sys.argv[0]] + args.split(" ")

    yield wrapper
    sys.argv = saved_argv


@pytest.fixture
def default_config():
    return Configuration()


@pytest.fixture
def fresh_tree(default_config):
    return SearchTree(config=default_config, root_smiles=None)


@pytest.fixture
def generate_root(default_config):
    def wrapper(smiles):
        return Node.create_root(smiles, tree=None, config=default_config)

    return wrapper


@pytest.fixture
def mock_policy(mocker, simple_actions):
    mocked_get_action = mocker.patch("aizynthfinder.mcts.policy.Policy.get_actions")

    def wrapper(mol):
        mocked_get_action.return_value = simple_actions(mol)
        return mocked_get_action.return_value

    return wrapper


@pytest.fixture
def mock_policy_model(mocker):
    class MockedKerasModel(mocker.MagicMock):
        @property
        def input(self):
            pass

        def predict(self, *_):
            pass

    mocker.patch.object(
        MockedKerasModel, "input", mocker.PropertyMock(return_value=np.zeros((3, 3)))
    )
    mocker.patch.object(
        MockedKerasModel,
        "predict",
        mocker.MagicMock(return_value=np.array([[0.2, 0.7, 0.1]])),
    )

    return mocker.patch(
        "aizynthfinder.utils.keras_utils.load_model", return_value=MockedKerasModel
    )


@pytest.fixture
def mock_stock(mocker):
    # This mock has to be applied after we expanded the root, due to the caching of the State methods
    def wrapper(availability_list):
        mocked_in_stock = mocker.patch("aizynthfinder.mcts.stock.Stock.__contains__")
        mocked_in_stock.side_effect = availability_list

    return wrapper


@pytest.fixture
def policy(default_config):
    policy = Policy(default_config)
    return policy


@pytest.fixture
def set_default_prior(default_config):
    default_config.use_prior = False

    def wrapper(prior):
        default_config.default_prior = prior

    yield wrapper
    default_config.use_prior = True


@pytest.fixture
def setup_search(fresh_tree, generate_root):
    tree = fresh_tree
    tree.initialize()

    def wrapper(smiles):
        root = generate_root(smiles)
        tree.add_root(root)
        return tree, root

    return wrapper


@pytest.fixture
def simple_actions():
    # These templated reactions are taken from the full USPTO data set
    # action1 and action2 can be applied, action3 cannot
    def wrapper(mol):
        actions = {
            "CCCCOc1ccc(CC(=O)N(C)O)cc1": [
                "([#8:4]-[N;H0;D3;+0:5](-[C;D1;H3:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])"
                ">>(Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([#8:4]-[NH;D2;+0:5]-[C;D1;H3:6])",
                "([C:2]-[CH2;D2;+0:1]-[O;H0;D2;+0:3]-[c:4])>>(Br-[CH2;D2;+0:1]-[C:2]).([OH;D1;+0:3]-[c:4])",
                "([C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])>>"
                "(O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([C:4]-[NH;D2;+0:5]-[C:6])",
            ]
        }
        action1, action2, action3 = actions[mol.smiles]
        action_list = [
            Reaction(mol, action1),
            Reaction(mol, action2),
            Reaction(mol, action3),
        ]
        prior_list = [0.7, 0.5, 0.3]
        return action_list, prior_list

    return wrapper


@pytest.fixture
def stock():
    stock = Stock()
    return stock


@pytest.fixture
def write_yaml(tmpdir):
    filename = str(tmpdir / "test.yaml")

    def wrapper(dict_):
        with open(filename, "w") as fileobj:
            yaml.dump(dict_, fileobj)
        return filename

    return wrapper
