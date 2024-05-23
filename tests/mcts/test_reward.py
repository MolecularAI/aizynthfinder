from aizynthfinder.analysis import TreeAnalysis
from aizynthfinder.context.scoring import (
    NumberOfReactionsScorer,
    StateScorer,
)
from aizynthfinder.search.mcts import MctsSearchTree


def test_reward_node(default_config, generate_root):
    config = default_config
    search_reward_scorer = repr(StateScorer(config))
    post_process_reward_scorer = repr(NumberOfReactionsScorer())

    config.search.algorithm_config["search_rewards"] = [search_reward_scorer]
    config.post_processing.route_scorers = [post_process_reward_scorer]

    node = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1", config)

    search_scorer = config.scorers[config.search.algorithm_config["search_rewards"][0]]
    route_scorer = config.scorers[config.post_processing.route_scorers[0]]

    assert round(search_scorer(node), 4) == 0.0491
    assert route_scorer(node) == 0


def test_default_postprocessing_reward(setup_aizynthfinder):
    """Test using default postprocessing.route_score"""
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    lookup = {root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0}}
    finder = setup_aizynthfinder(lookup, child1_smi)

    config = finder.config
    config.search.return_first = True

    search_reward_scorer = repr(NumberOfReactionsScorer())
    state_scorer = repr(StateScorer(config))

    config.search.algorithm_config["search_rewards"] = [search_reward_scorer]
    finder.config = config

    assert len(finder.config.post_processing.route_scorers) == 0

    finder.tree_search()
    tree_analysis_search = TreeAnalysis(
        finder.tree, scorer=config.scorers[search_reward_scorer]
    )
    tree_analysis_pp = TreeAnalysis(finder.tree, scorer=config.scorers[state_scorer])

    finder.build_routes()
    assert finder.tree.reward_scorer_name == search_reward_scorer

    top_score_tree_analysis = tree_analysis_search.tree_statistics()["top_score"]
    top_score_finder = finder.tree.compute_reward(tree_analysis_search.best())

    assert top_score_finder == top_score_tree_analysis

    top_score_tree_analysis = tree_analysis_pp.tree_statistics()["top_score"]
    top_score_finder = finder.analysis.tree_statistics()["top_score"]

    # Finder used the search_reward_scorer and not state_scorer
    assert top_score_finder != top_score_tree_analysis


def test_custom_reward(setup_aizynthfinder):
    """Test using different custom reward functions for MCTS and route building."""

    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    lookup = {root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0}}
    finder = setup_aizynthfinder(lookup, child1_smi)

    # Test first with return_first and multiple route scores
    config = finder.config
    config.search.return_first = True

    search_reward_scorer = repr(StateScorer(config))
    post_process_reward_scorer = repr(NumberOfReactionsScorer())

    config.search.algorithm_config["search_rewards"] = [search_reward_scorer]
    config.post_processing.route_scorers = [post_process_reward_scorer]
    finder.config = config

    assert finder.config.post_processing.route_scorers == [post_process_reward_scorer]

    finder.tree_search()
    tree_analysis_search = TreeAnalysis(
        finder.tree, scorer=config.scorers[search_reward_scorer]
    )
    tree_analysis_pp = TreeAnalysis(
        finder.tree, scorer=config.scorers[post_process_reward_scorer]
    )

    finder.build_routes()

    assert finder.config.post_processing.route_scorers == [post_process_reward_scorer]
    assert finder.tree.reward_scorer_name == search_reward_scorer

    top_score_tree_analysis = tree_analysis_search.tree_statistics()["top_score"]
    top_score_finder = finder.tree.compute_reward(tree_analysis_search.best())

    assert top_score_finder == top_score_tree_analysis

    top_score_tree_analysis = tree_analysis_pp.tree_statistics()["top_score"]
    top_score_finder = finder.analysis.tree_statistics()["top_score"]

    assert top_score_finder == top_score_tree_analysis


def test_reward_node_backward_compatibility(default_config):
    reward_scorer = repr(NumberOfReactionsScorer())
    default_config.search.algorithm_config["search_reward"] = reward_scorer

    tree = MctsSearchTree(config=default_config, root_smiles=None)

    assert tree.reward_scorer_name == reward_scorer
