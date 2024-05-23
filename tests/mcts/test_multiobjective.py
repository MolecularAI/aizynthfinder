import pytest

from aizynthfinder.search.mcts.node import ParetoMctsNode
from aizynthfinder.search.mcts import MctsSearchTree


@pytest.fixture
def generate_root(default_config):
    def wrapper(smiles, config=None):
        return ParetoMctsNode.create_root(
            smiles, tree=None, config=config or default_config
        )

    return wrapper


@pytest.fixture
def setup_mcts_search(
    default_config, get_one_step_expansion, setup_policies, generate_root
):
    default_config.search.algorithm_config["search_rewards"] = [
        "number of reactions",
        "number of pre-cursors in stock",
    ]
    expansion_strategy, filter_strategy = setup_policies(get_one_step_expansion)
    root_smiles = list(expansion_strategy.lookup.keys())[0]
    return (
        generate_root(root_smiles),
        expansion_strategy,
        filter_strategy,
    )


def test_expand_root_node(setup_mcts_search):
    root, _, _ = setup_mcts_search

    root.expand()

    view = root.children_view()
    assert len(view["actions"]) == 3
    assert view["priors"] == [[0.7, 0.7], [0.5, 0.5], [0.3, 0.3]]
    assert view["values"] == [[0.7, 0.7], [0.5, 0.5], [0.3, 0.3]]
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_expand_root_with_default_priors(setup_mcts_search, set_default_prior):
    root, _, _ = setup_mcts_search
    set_default_prior(0.01)

    root.expand()

    view = root.children_view()
    assert len(view["actions"]) == 3
    assert view["priors"] == [[0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]
    assert view["values"] == [[0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_backpropagate(setup_mcts_search):
    root, _, _ = setup_mcts_search
    root.expand()
    child = root.promising_child()
    view_prior = root.children_view()

    root.backpropagate(child, [1.5, 2.0])

    view_post = root.children_view()
    assert view_post["visitations"][0] == view_prior["visitations"][0] + 1
    assert view_prior["visitations"][1:] == view_post["visitations"][1:]
    assert view_prior["values"] == view_post["values"]
    assert view_post["rewards_cum"][0][0] == view_prior["rewards_cum"][0][0] + 1.5
    assert view_post["rewards_cum"][0][1] == view_prior["rewards_cum"][0][1] + 2.0
    assert view_prior["rewards_cum"][1:] == view_post["rewards_cum"][1:]


def test_setup_weighted_sum_tree(default_config):
    default_config.search.algorithm_config["search_rewards"] = [
        "number of reactions",
        "number of pre-cursors in stock",
    ]
    default_config.search.algorithm_config["search_rewards_weights"] = [0.5, 0.1]
    root_smiles = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    tree = MctsSearchTree(config=default_config, root_smiles=root_smiles)

    assert tree.mode == "weighted-sum"
    assert len(tree.reward_scorer.selection) == 2
    assert tree.compute_reward(tree.root) == 0.0


def test_setup_mo_tree(default_config):
    default_config.search.algorithm_config["search_rewards"] = [
        "number of reactions",
        "number of pre-cursors in stock",
    ]
    root_smiles = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    tree = MctsSearchTree(config=default_config, root_smiles=root_smiles)

    assert tree.mode == "multi-objective"
    assert len(tree.reward_scorer.selection) == 2
    assert tree.compute_reward(tree.root) == [0.0, 0.0]
