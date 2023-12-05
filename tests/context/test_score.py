import numpy as np
import pytest
from aizynthfinder.chem import Molecule, UniqueMolecule
from aizynthfinder.context.scoring import (
    AverageTemplateOccurrenceScorer,
    CombinedScorer,
    FractionInStockScorer,
    MaxTransformScorerer,
    NumberOfPrecursorsInStockScorer,
    NumberOfPrecursorsScorer,
    NumberOfReactionsScorer,
    PriceSumScorer,
    ReactionClassMembershipScorer,
    RouteCostScorer,
    ScorerCollection,
    ScorerException,
    StockAvailabilityScorer,
    StateScorer,
)
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsSearchTree


def test_state_scorer_node(default_config, setup_linear_mcts):
    _, node = setup_linear_mcts()
    scorer = StateScorer(default_config)

    assert repr(scorer) == "state score"
    assert round(scorer(node), 4) == 0.994


def test_state_scorer_nodes(setup_linear_mcts, setup_branched_mcts, default_config):
    _, node1 = setup_linear_mcts()
    _, node2 = setup_branched_mcts()
    scorer = StateScorer(default_config)

    scores = scorer([node1, node2])

    assert repr(scorer) == "state score"
    assert round(scores[0], 4) == 0.994
    assert round(scores[1], 4) == 0.9866


def test_state_scorer_tree(default_config, setup_linear_reaction_tree):
    tree = setup_linear_reaction_tree()
    scorer = StateScorer(default_config)

    assert round(scorer(tree), 4) == 0.994


def test_state_scorer_trees(default_config, setup_linear_reaction_tree):
    rt = setup_linear_reaction_tree()
    scorer = StateScorer(default_config)

    scores = scorer([rt, rt])

    assert round(scores[0], 4) == 0.994
    assert round(scores[1], 4) == 0.994


def test_sort(default_config, setup_linear_mcts, setup_branched_mcts):
    _, node1 = setup_linear_mcts()
    _, node2 = setup_branched_mcts()
    scorer = StateScorer(default_config)

    sorted_nodes, scores, _ = scorer.sort([node2, node1])

    assert [np.round(score, 4) for score in scores] == [0.994, 0.9866]
    assert sorted_nodes == [node1, node2]


def test_template_occurrence_scorer_no_metadata(setup_linear_mcts, default_config):
    _, node1 = setup_linear_mcts()
    scorer = AverageTemplateOccurrenceScorer(config=default_config)

    assert scorer(node1) == 0


def test_template_occurrence_scorer(setup_linear_mcts, default_config):
    search_tree, _ = setup_linear_mcts()
    nodes = list(search_tree.graph())
    nodes[0][nodes[1]]["action"].metadata["library_occurrence"] = 5
    nodes[1][nodes[2]]["action"].metadata["library_occurence"] = 10
    scorer = AverageTemplateOccurrenceScorer(default_config)

    assert scorer(nodes[0]) == 0
    assert scorer(nodes[1]) == 5
    assert scorer(nodes[2]) == 7.5


def test_template_occurrence_scorer_minmax_scaled(setup_linear_mcts, default_config):
    search_tree, _ = setup_linear_mcts()
    nodes = list(search_tree.graph())
    nodes[0][nodes[1]]["action"].metadata["library_occurrence"] = 5
    nodes[1][nodes[2]]["action"].metadata["library_occurence"] = 10

    scaling_params = {"name": "min_max", "min_val": 1, "max_val": 10, "reverse": False}
    scorer = AverageTemplateOccurrenceScorer(default_config, scaling_params)
    scores = scorer(nodes[:3])
    assert [round(score, 1) for score in scores] == [0, 0.4, 0.7]


def test_template_occurrence_scorer_squash_scaled(setup_linear_mcts, default_config):
    search_tree, _ = setup_linear_mcts()
    nodes = list(search_tree.graph())
    nodes[0][nodes[1]]["action"].metadata["library_occurrence"] = 5
    nodes[1][nodes[2]]["action"].metadata["library_occurence"] = 10

    scaling_params = {"name": "squash", "slope": 1, "xoffset": 1, "yoffset": 0}
    scorer = AverageTemplateOccurrenceScorer(default_config, scaling_params)
    scores = scorer(nodes[:3])
    assert [round(score, 4) for score in scores] == [0.2689, 0.9820, 0.9985]


def test_template_occurrence_scorer_tree(setup_linear_reaction_tree, default_config):
    tree = setup_linear_reaction_tree()
    scorer = AverageTemplateOccurrenceScorer(default_config)

    assert scorer(tree) == 0


def test_template_occurrence_scorer_tree_one_node(default_config):
    rt = ReactionTree()
    rt.root = Molecule(smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    rt.graph.add_node(rt.root)
    scorer = AverageTemplateOccurrenceScorer(default_config)

    assert scorer(rt) == 0.0


def test_scorers_one_mcts_node(default_config):
    tree = MctsSearchTree(default_config, root_smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    node = tree.root

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-3) == 0.0497
    assert FractionInStockScorer(default_config)(node) == 0
    assert MaxTransformScorerer(default_config)(node) == 0
    assert NumberOfReactionsScorer(default_config)(node) == 0
    assert NumberOfPrecursorsScorer(default_config)(node) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 0
    assert PriceSumScorer(default_config)(node) == 10
    assert RouteCostScorer(default_config)(node) == 10


def test_scoring_branched_mcts_tree(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-4) == 0.9866
    assert FractionInStockScorer(default_config)(node) == 1
    assert MaxTransformScorerer(default_config)(node) == 3
    assert NumberOfReactionsScorer()(node) == 4
    assert NumberOfPrecursorsScorer(default_config)(node) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 5
    assert PriceSumScorer(default_config)(node) == 5
    cost_score = RouteCostScorer(default_config)(node)
    assert pytest.approx(cost_score, abs=1e-4) == 13.6563


def test_scoring_branched_mcts_tree_not_in_stock(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts("O")

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-4) == 0.7966
    assert FractionInStockScorer(default_config)(node) == 4 / 5
    assert MaxTransformScorerer(default_config)(node) == 3
    assert NumberOfReactionsScorer()(node) == 4
    assert NumberOfPrecursorsScorer(default_config)(node) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 4
    assert PriceSumScorer(default_config)(node) == 14
    cost_score = RouteCostScorer(default_config)(node)
    assert pytest.approx(cost_score, abs=1e-4) == 31.2344


def test_scorers_tree_one_node_route(default_config):
    tree = ReactionTree()
    tree.root = UniqueMolecule(smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    tree.graph.add_node(tree.root)

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-3) == 0.0497
    assert FractionInStockScorer(default_config)(tree) == 0
    assert MaxTransformScorerer(default_config)(tree) == -1
    assert NumberOfReactionsScorer(default_config)(tree) == 0
    assert NumberOfPrecursorsScorer(default_config)(tree) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 0
    assert PriceSumScorer(default_config)(tree) == 10
    assert RouteCostScorer(default_config)(tree) == 10


def test_scoring_branched_route(default_config, setup_branched_reaction_tree):
    tree = setup_branched_reaction_tree()

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-4) == 0.9866
    assert FractionInStockScorer(default_config)(tree) == 1
    assert MaxTransformScorerer(default_config)(tree) == 3
    assert NumberOfReactionsScorer()(tree) == 4
    assert NumberOfPrecursorsScorer(default_config)(tree) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 5
    assert PriceSumScorer(default_config)(tree) == 5
    cost_score = RouteCostScorer(default_config)(tree)
    assert pytest.approx(cost_score, abs=1e-4) == 13.6563


def test_scoring_branched_route_not_in_stock(
    default_config, setup_branched_reaction_tree
):
    tree = setup_branched_reaction_tree("O")

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-4) == 0.7966
    assert FractionInStockScorer(default_config)(tree) == 4 / 5
    assert MaxTransformScorerer(default_config)(tree) == 3
    assert NumberOfReactionsScorer()(tree) == 4
    assert NumberOfPrecursorsScorer(default_config)(tree) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 4
    assert PriceSumScorer(default_config)(tree) == 14
    cost_score = RouteCostScorer(default_config)(tree)
    assert pytest.approx(cost_score, abs=1e-4) == 31.2344


def test_reaction_class_scorer_tree(default_config, load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))
    scorer = ReactionClassMembershipScorer(default_config, reaction_class_set=["abc"])

    # one within the set, one outside the set
    assert scorer(tree) == 0.1


def test_stock_availability_scorer_tree(
    default_config, load_reaction_tree, setup_stock
):
    setup_stock(
        default_config,
        "NC1CCCC(C2C=CC=C2)C1",
        "c1ccccc1",
    )
    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))
    scorer = StockAvailabilityScorer(default_config, source_score={"stock": 0.5})

    # 0.5*0.5*0.1 two in stock, one not in stock
    assert scorer(tree) == 0.025


def test_create_scorer_collection(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()

    assert len(collection) == 5

    assert "state score" in collection.names()
    assert "number of reactions" in collection.names()

    assert isinstance(collection["state score"], StateScorer)

    with pytest.raises(KeyError):
        collection["dummy"]


def test_delete_scorer_to_collection(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()

    del collection["state score"]

    assert "state score" not in collection.names()


def test_add_scorer_to_collection(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()
    del collection["state score"]

    collection.load(StateScorer(default_config))

    assert "state score" in collection.names()


def test_add_scorer_to_collection_no_scorer(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()

    with pytest.raises(ScorerException):
        collection.load(Molecule(smiles="CCC"))


def test_load_scorer_to_collection_only_class(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()
    del collection["state score"]

    collection.load_from_config(**{"StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_full_package(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()
    del collection["state score"]

    collection.load_from_config(**{"aizynthfinder.context.scoring.StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_failures(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()

    with pytest.raises(ScorerException, match=".*load module.*"):
        collection.load_from_config(**{"mypackage.scoring.StateScorer": {}})

    with pytest.raises(ScorerException, match=".*class.*"):
        collection.load_from_config(**{"aizynthfinder.context.scoring.NoScorer": {}})


def test_subset_scorer_collection(default_config):
    collection = ScorerCollection(default_config)
    names = ["state score", "number of pre-cursors"]
    new_collection = collection.make_subset(names)

    assert len(new_collection) == 2
    assert new_collection.names() == names


def test_subset_scorer_collection_not_exists(default_config):
    collection = ScorerCollection(default_config)
    names = ["state score", "dummy"]

    with pytest.raises(ScorerException, match="dummy"):
        collection.make_subset(names)


def test_score_vector(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)
    collection.select_all()

    scores = collection.score_vector(node)
    assert pytest.approx(scores, abs=1e-3) == [0.9866, 4, 5, 5, 0]


def test_score_vector_no_selection(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)

    scores = collection.score_vector(node)
    assert scores == []


def test_weighted_sum_score(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)
    collection.select_all()

    score = collection.weighted_score(node, weights=[0.0, 0.5, 1.0, 1.0, 0.0])
    assert score == 12

    with pytest.raises(ScorerException):
        collection.weighted_score(node, weights=[])

    with pytest.raises(ScorerException):
        collection.weighted_score(node, weights=[1.0])


def test_weighted_sum_score_no_selection(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)

    with pytest.raises(ScorerException):
        collection.weighted_score(node, weights=[0.0, 0.5, 1.0, 1.0, 0.0])


def test_combined_scorer_node_default_weights(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()

    combined_scorer = CombinedScorer(
        default_config, ["state score", "number of reactions"]
    )
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + number of reactions"
    assert round(score, 4) == 2.4933


def test_combined_scorer_node(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()

    combined_scorer = CombinedScorer(
        default_config,
        ["state score", "number of reactions"],
        [0.75, 0.25],
    )
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + number of reactions"
    assert round(score, 4) == 1.7399


def test_combined_scorer_reaction_tree(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    reaction_tree = node.to_reaction_tree()

    combined_scorer = CombinedScorer(
        default_config,
        ["state score", "number of reactions"],
        [0.75, 0.25],
    )
    score = combined_scorer(reaction_tree)

    assert repr(combined_scorer) == "state score + number of reactions"
    assert round(score, 4) == 1.7399
