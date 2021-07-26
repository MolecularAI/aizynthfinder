import pytest

from aizynthfinder.search.mcts import MctsNode, MctsSearchTree


@pytest.fixture
def generate_root(default_config):
    def wrapper(smiles, config=None):
        return MctsNode.create_root(smiles, tree=None, config=config or default_config)

    return wrapper


@pytest.fixture
def fresh_tree(default_config):
    return MctsSearchTree(config=default_config, root_smiles=None)


@pytest.fixture
def set_default_prior(default_config):
    default_config.use_prior = False

    def wrapper(prior):
        default_config.default_prior = prior

    yield wrapper
    default_config.use_prior = True


@pytest.fixture
def setup_mcts_search(get_one_step_expansion, setup_policies, generate_root):
    expansion_strategy, filter_strategy = setup_policies(get_one_step_expansion)
    root_smiles = list(expansion_strategy.lookup.keys())[0]
    return (
        generate_root(root_smiles),
        expansion_strategy,
        filter_strategy,
    )


@pytest.fixture
def setup_complete_mcts_tree(default_config, setup_policies, setup_stock):
    root_smiles = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    tree = MctsSearchTree(config=default_config, root_smiles=root_smiles)
    lookup = {
        root_smiles: {
            "smiles": "CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O",
            "prior": 1.0,
        },
        "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F": {
            "smiles": "N#Cc1cccc(N)c1F.O=C(Cl)c1ccc(F)cc1",
            "prior": 1.0,
        },
    }
    setup_policies(lookup)

    setup_stock(
        default_config, "CN1CCC(Cl)CC1", "O", "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"
    )

    node1 = tree.root
    node1.expand()

    node2 = node1.promising_child()
    node2.expand()

    node3 = node2.promising_child()
    node3.expand()

    return tree, [node1, node2, node3]
