import pytest
import numpy as np


from aizynthfinder.search.retrostar.search_tree import SearchTree
from aizynthfinder.search.retrostar.nodes import MoleculeNode
from aizynthfinder.aizynthfinder import AiZynthFinder


@pytest.fixture
def setup_aizynthfinder(setup_policies, setup_stock):
    def wrapper(expansions, stock):
        finder = AiZynthFinder()
        root_smi = list(expansions.keys())[0]
        setup_policies(expansions, config=finder.config)
        setup_stock(finder.config, *stock)
        finder.target_smiles = root_smi
        finder.config.search_algorithm = (
            "aizynthfinder.search.retrostar.search_tree.SearchTree"
        )
        return finder

    return wrapper


@pytest.fixture
def setup_mocked_model(mocker):
    biases = [np.zeros(10), np.zeros(1)]
    weights = [np.ones([10, 10]), np.ones([10, 1])]

    mocker.patch("builtins.open")
    mocked_pickle_load = mocker.patch("aizynthfinder.search.retrostar.cost.pickle.load")
    mocked_pickle_load.return_value = weights, biases


@pytest.fixture
def setup_search_tree(default_config, setup_policies, setup_stock):
    root_smiles = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    tree = SearchTree(config=default_config, root_smiles=root_smiles)
    lookup = {
        root_smiles: {
            "smiles": "CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O",
            "prior": 1.0,
        }
    }
    setup_policies(lookup)

    setup_stock(
        default_config, "CN1CCC(Cl)CC1", "O", "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"
    )
    return tree


@pytest.fixture
def setup_star_root(default_config):
    def wrapper(smiles):
        return MoleculeNode.create_root(smiles, config=default_config)

    return wrapper
