from networkx import get_node_attributes, set_node_attributes
import numpy as np
import pytest


from aizynthfinder.chem import Molecule, UniqueMolecule
from aizynthfinder.context.config import Configuration
from aizynthfinder.context.scoring import (
    AverageTemplateOccurrenceScorer,
    BrokenBondsScorer,
    CombinedScorer,
    DeepSetScorer,
    DeltaSyntheticComplexityScorer,
    FractionInSourceStockScorer,
    FractionInStockScorer,
    MaxTransformScorer,
    FractionOfIntermediatesInStockScorer,
    NumberOfPrecursorsInStockScorer,
    NumberOfPrecursorsScorer,
    NumberOfReactionsScorer,
    PriceSumScorer,
    ReactionClassMembershipScorer,
    ReactionClassRankScorer,
    RouteCostScorer,
    RouteSimilarityScorer,
    ScorerCollection,
    ScorerException,
    StockAvailabilityScorer,
    StateScorer,
    SUPPORT_DISTANCES,
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
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(node) == 0
    )
    assert MaxTransformScorer(default_config)(node) == 0
    assert NumberOfReactionsScorer(default_config)(node) == 0
    assert NumberOfPrecursorsScorer(default_config)(node) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 0
    assert PriceSumScorer(default_config)(node) == 10
    assert RouteCostScorer(default_config)(node) == 10


def test_scoring_branched_mcts_tree(default_config, setup_branched_mcts, setup_stock):
    _, node = setup_branched_mcts()

    config2 = Configuration()
    setup_stock(config2, "NC1CCCC(C2C=CC=C2)C1", "OOc1ccccc1")
    default_config.stock.load(config2.stock["stock"], "stock2")
    default_config.stock.select(["stock"])

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-4) == 0.9866
    assert FractionInStockScorer(default_config)(node) == 1
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(node) == 1
    )
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["not-a-stock"])(node)
        == 0
    )

    assert (
        round(
            FractionOfIntermediatesInStockScorer(default_config, stock_name="stock2")(
                node
            ),
            4,
        )
        == 0.6667
    )
    assert MaxTransformScorer(default_config)(node) == 3
    assert NumberOfReactionsScorer()(node) == 4
    assert NumberOfPrecursorsScorer(default_config)(node) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 5
    assert PriceSumScorer(default_config)(node) == 5
    cost_score = RouteCostScorer(default_config)(node)
    assert pytest.approx(cost_score, abs=1e-4) == 13.6563


def test_scoring_branched_mcts_tree_root(default_config, setup_stock):
    # Test a route that is solved but has no reactions
    setup_stock(None, "CCCCOc1ccc(CC(=O)N(C)O)cc1")
    tree = MctsSearchTree(default_config, root_smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    node = tree.root

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-4) == 0.9991
    assert FractionInStockScorer(default_config)(node) == 1.0
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(node)
        == 1.0
    )
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["not-a-stock"])(node)
        == 0
    )
    assert (
        FractionOfIntermediatesInStockScorer(default_config, stock_name="stock")(node)
        == 1.0
    )
    assert MaxTransformScorer(default_config)(node) == 0
    assert NumberOfReactionsScorer()(node) == 0
    assert NumberOfPrecursorsScorer(default_config)(node) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 1
    assert PriceSumScorer(default_config)(node) == 1.0
    assert RouteCostScorer(default_config)(node) == 1.0
    assert StockAvailabilityScorer(default_config, {"stock": 1})(node) == 1.0
    assert ReactionClassMembershipScorer(default_config, ["abc"])(node) == 1.0


def test_scoring_branched_mcts_tree_not_in_stock(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts("O")

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-4) == 0.7966
    assert FractionInStockScorer(default_config)(node) == 4 / 5
    assert MaxTransformScorer(default_config)(node) == 3
    assert NumberOfReactionsScorer()(node) == 4
    assert NumberOfPrecursorsScorer(default_config)(node) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 4
    assert PriceSumScorer(default_config)(node) == 14
    cost_score = RouteCostScorer(default_config)(node)
    assert pytest.approx(cost_score, abs=1e-4) == 31.2344


def test_scorers_tree_one_node_route(default_config, setup_stock, shared_datadir):
    setup_stock(None, "CCCCOc1ccc(CC(=O)N(C)O)cc1")
    tree = ReactionTree()
    tree.root = UniqueMolecule(smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    tree.graph.add_node(tree.root)
    tree.graph.nodes[tree.root]["in_stock"] = True

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-3) == 0.9991
    assert FractionInStockScorer(default_config)(tree) == 1
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(tree)
        == 1.0
    )
    assert (
        FractionOfIntermediatesInStockScorer(default_config, stock_name="stock")(tree)
        == 1.0
    )
    assert MaxTransformScorer(default_config)(tree) == -1
    assert NumberOfReactionsScorer(default_config)(tree) == 0
    assert NumberOfPrecursorsScorer(default_config)(tree) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 1
    assert PriceSumScorer(default_config)(tree) == 1
    assert RouteCostScorer(default_config)(tree) == 1
    assert StockAvailabilityScorer(default_config, {"stock": 1})(tree) == 1.0
    assert ReactionClassMembershipScorer(default_config, ["abc"])(tree) == 1.0


def test_scoring_branched_route(default_config, setup_branched_reaction_tree):
    tree = setup_branched_reaction_tree()

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-4) == 0.9866
    assert FractionInStockScorer(default_config)(tree) == 1
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(tree) == 1
    )
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["not-a-stock"])(tree)
        == 0
    )
    assert MaxTransformScorer(default_config)(tree) == 3
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
    assert (
        FractionInSourceStockScorer(default_config, source_stocks=["stock"])(tree)
        == 4 / 5
    )
    assert FractionInStockScorer(default_config)(tree) == 4 / 5
    assert MaxTransformScorer(default_config)(tree) == 3
    assert NumberOfReactionsScorer()(tree) == 4
    assert NumberOfPrecursorsScorer(default_config)(tree) == 5
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 4
    assert PriceSumScorer(default_config)(tree) == 14
    cost_score = RouteCostScorer(default_config)(tree)
    assert pytest.approx(cost_score, abs=1e-4) == 31.2344



def test_reaction_class_scorer_tree(shared_datadir, default_config, load_reaction_tree):
    tree = ReactionTree.from_dict(
        load_reaction_tree("linear_route_w_metadata.json", remove_metadata=False)
    )
    scorer = ReactionClassMembershipScorer(default_config, reaction_class_set=["abc"])

    # one within the set, one outside the set
    assert round(scorer(tree), 4) == 0.1
    scorer2 = ReactionClassMembershipScorer(
        default_config, reaction_class_set=str(shared_datadir / "reaction_classes.txt")
    )
    assert round(scorer2(tree), 4) == 0.1

    # both reaction classes in the set
    scorer3 = ReactionClassMembershipScorer(
        default_config, reaction_class_set=["abc", "xyz"]
    )
    assert round(scorer3(tree), 4) == 1.0

    # empty list => no reaction classes in the set
    scorer4 = ReactionClassMembershipScorer(default_config, reaction_class_set=[])
    assert round(scorer4(tree), 4) == 0.01


def test_reaction_scorer_w_power_scaling(default_config, load_reaction_tree):
    n_reaction_scorer = NumberOfReactionsScorer(default_config)
    scorer = NumberOfReactionsScorer(default_config, scaler_params={"name": "power"})

    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))

    assert n_reaction_scorer(tree) == 2
    assert round(scorer(tree), 4) == 0.9604

    tree = ReactionTree.from_dict(load_reaction_tree("branched_route.json"))

    assert n_reaction_scorer(tree) == 4
    assert round(scorer(tree), 4) == 0.9224


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


def test_stock_availability_scorer_tree_other_source_score(
    default_config, load_reaction_tree, setup_stock
):
    # Setting up two stocks for this test
    setup_stock(
        default_config,
        "NC1CCCC(C2C=CC=C2)C1",
    )
    config2 = Configuration()
    setup_stock(
        config2,
        "c1ccccc1",
    )
    default_config.stock.load(config2.stock["stock"], "stock2")
    default_config.stock.select(["stock", "stock2"])
    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))
    # Only specifying one stock explicitly
    scorer = StockAvailabilityScorer(default_config, source_score={"stock": 0.5})

    # 0.5*0.1*0.1 one in stock, one in unspecified stock, one not in stock
    assert pytest.approx(scorer(tree)) == 0.005

    scorer.other_source_score = 0.2

    # 0.5*0.2*0.1 one in stock, one in unspecified stock, one not in stock
    assert pytest.approx(scorer(tree)) == 0.01


def test_fraction_in_source_stock_score(
    default_config, load_reaction_tree, setup_stock
):
    # Setting up two stocks for this test
    setup_stock(
        default_config,
        "NC1CCCC(C2C=CC=C2)C1",
    )
    config2 = Configuration()
    setup_stock(
        config2,
        "c1ccccc1",
    )
    default_config.stock.load(config2.stock["stock"], "stock2")
    default_config.stock.select(["stock", "stock2"])
    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))

    # Not specifying source stock
    with pytest.raises(
        TypeError, match="missing 1 required positional argument: 'source_stocks'"
    ):
        FractionInSourceStockScorer(default_config)

    # Only specifying one stock explicitly
    scorer = FractionInSourceStockScorer(default_config, source_stocks=["stock2"])

    # (0+1+0)/3 one in stock, one in stock2, one not in stock
    assert round(scorer(tree), 4) == 0.3333

    # Specifying both stocks explicitly
    scorer = FractionInSourceStockScorer(
        default_config, source_stocks=["stock", "stock2"]
    )

    # (1+1+0)/3 one in stock, one in stock2, one not in stock
    assert round(scorer(tree), 4) == 0.6667


def test_fraction_intermediates_in_stock_score(
    default_config, load_reaction_tree, setup_stock
):
    # Setting up two stocks for this test
    setup_stock(
        default_config,
        "NC1CCCC(C2C=CC=C2)C1",
        "c1ccccc1",
    )
    config2 = Configuration()
    setup_stock(config2, "OOc1ccc(-c2ccccc2)cc1")
    default_config.stock.load(config2.stock["stock"], "stock2")
    default_config.stock.select(["stock"])
    tree = ReactionTree.from_dict(load_reaction_tree("linear_route_w_metadata.json"))

    # Not specifying source stock
    with pytest.raises(
        TypeError, match="missing 1 required positional argument: 'stock_name'"
    ):
        FractionOfIntermediatesInStockScorer(default_config)

    # Only specifying one stock explicitly
    scorer = FractionOfIntermediatesInStockScorer(default_config, stock_name="stock2")

    # one intermediate in stock2 (the only intermediate node which is not in stock)
    assert round(scorer(tree), 4) == 1

    # Specifying stock the same stock as the one used in search
    scorer = FractionOfIntermediatesInStockScorer(default_config, stock_name="stock")

    # no intermediates are found in stock since the same stock was used in the search
    assert round(scorer(tree), 4) == 0


def test_create_scorer_collection(default_config):
    collection = ScorerCollection(default_config)
    collection.create_default_scorers()

    assert len(collection) == 4

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
    assert pytest.approx(scores, abs=1e-3) == [0.9866, 4, 5, 5]


def test_score_vector_no_selection(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)

    scores = collection.score_vector(node)
    assert scores == []


def test_weighted_sum_score(default_config, setup_branched_mcts):
    _, node = setup_branched_mcts()
    collection = ScorerCollection(default_config)
    collection.select_all()

    score = collection.weighted_score(node, weights=[0.0, 0.5, 1.0, 1.0])
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


def test_broken_bonds_scorer_node(default_config, setup_mcts_broken_bonds):
    _, node = setup_mcts_broken_bonds()

    default_config.search.break_bonds = [(1, 2)]
    scorer = BrokenBondsScorer(default_config)(node)

    assert scorer == 0.5


def test_broken_bonds_scorer_node_unbroken(default_config, setup_mcts_broken_bonds):
    _, node = setup_mcts_broken_bonds(broken=False)

    default_config.search.break_bonds = [(1, 2)]
    scorer = BrokenBondsScorer(default_config)(node)

    assert scorer == 0.0


def test_broken_bonds_or_operator_scorer_node_unbroken(
    default_config, setup_mcts_broken_bonds
):
    _, node = setup_mcts_broken_bonds(broken=False)

    default_config.search.break_bonds = [(1, 2)]
    default_config.search.break_bonds_operator = "or"
    scorer = BrokenBondsScorer(default_config)(node)

    assert scorer == 0.0


def test_broken_bonds_scorer_reaction_tree(default_config, setup_mcts_broken_bonds):
    _, node = setup_mcts_broken_bonds()
    reaction_tree = node.to_reaction_tree()

    default_config.search.break_bonds = [(1, 2)]
    scorer = BrokenBondsScorer(default_config)(reaction_tree)

    assert round(scorer, 4) == 0.5


def test_combined_scorer_node_default_weights(default_config, setup_mcts_broken_bonds):
    default_config.search.break_bonds = [(1, 2)]
    _, node = setup_mcts_broken_bonds(config=default_config)

    combined_scorer = CombinedScorer(default_config, ["state score", "broken bonds"])
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + broken bonds"
    assert round(score, 4) == 0.747


def test_combined_scorer_node_default_weights_geometric(
    default_config, setup_mcts_broken_bonds
):
    default_config.search.break_bonds = [(1, 2)]
    _, node = setup_mcts_broken_bonds(config=default_config)

    combined_scorer = CombinedScorer(
        default_config,
        ["state score", "broken bonds"],
        combine_strategy="mean-geometric",
    )
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + broken bonds (mean-geometric)"
    assert round(score, 4) == 0.705


def test_combined_scorer_node_default_weights_product(
    default_config, setup_mcts_broken_bonds
):
    default_config.search.break_bonds = [(1, 2)]
    _, node = setup_mcts_broken_bonds(config=default_config)

    combined_scorer = CombinedScorer(
        default_config, ["state score", "broken bonds"], combine_strategy="product"
    )
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + broken bonds (product)"
    assert round(score, 4) == 0.497


def test_combined_scorer_node(default_config, setup_mcts_broken_bonds):
    default_config.search.break_bonds = [(1, 2)]
    _, node = setup_mcts_broken_bonds(config=default_config)

    combined_scorer = CombinedScorer(
        default_config,
        ["state score", "broken bonds", "number of reactions"],
        [0.48, 0.5, 0.02],
    )
    score = combined_scorer(node)

    assert repr(combined_scorer) == "state score + broken bonds + number of reactions"
    assert round(score, 4) == 0.7671


def test_combined_scorer_reaction_tree(default_config, setup_mcts_broken_bonds):
    default_config.search.break_bonds = [(1, 2)]
    _, node = setup_mcts_broken_bonds(config=default_config)
    reaction_tree = node.to_reaction_tree()

    combined_scorer = CombinedScorer(
        default_config,
        ["state score", "broken bonds", "number of reactions"],
        [0.48, 0.5, 0.02],
    )
    score = combined_scorer(reaction_tree)

    assert repr(combined_scorer) == "state score + broken bonds + number of reactions"
    assert round(score, 4) == 0.7671


@pytest.mark.xfail(
    condition=not SUPPORT_DISTANCES, reason="route_distance package not installed"
)
def test_route_similarity_no_ref_routes(
    default_config, mocker, setup_branched_reaction_tree
):
    tree = setup_branched_reaction_tree()
    mocker.patch("aizynthfinder.context.scoring.scorers.route_distances_calculator")

    scorer = RouteSimilarityScorer(default_config, "", "dummy")

    assert scorer(tree) == 1.0


@pytest.mark.xfail(
    condition=not SUPPORT_DISTANCES, reason="route_distance package not installed"
)
def test_route_similarity_self(default_config, mocker, setup_branched_reaction_tree):
    tree = setup_branched_reaction_tree()
    calc_patch = mocker.patch(
        "aizynthfinder.context.scoring.scorers.route_distances_calculator"
    )
    calc_patch.return_value.return_value = np.asarray([[0.0, 0.0], [0.0, 0.0]])

    scorer = RouteSimilarityScorer(default_config, "", "dummy")
    scorer.routes = [tree]
    scorer.n_routes = 1

    assert pytest.approx(scorer(tree), abs=1e-2) == 0.0


@pytest.mark.xfail(
    condition=not SUPPORT_DISTANCES, reason="route_distance package not installed"
)
def test_route_similarity_self_self(
    default_config, mocker, setup_branched_reaction_tree
):
    tree = setup_branched_reaction_tree()
    calc_patch = mocker.patch(
        "aizynthfinder.context.scoring.scorers.route_distances_calculator"
    )
    calc_patch.return_value.return_value = np.asarray(
        [[10.0, 0.0, 0.0], [0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
    )

    scorer = RouteSimilarityScorer(default_config, "", "dummy")
    scorer.routes = [tree, tree]
    scorer.n_routes = 2

    assert pytest.approx(scorer(tree), abs=1e-3) == 0.007

    scorer.similarity = True

    assert pytest.approx(scorer(tree), abs=1e-3) == 0.993

    scorer.agg_func = np.max

    assert pytest.approx(scorer(tree), abs=1e-2) == 0.5


def test_delta_complexity_scorer_tree(
    default_config, mocker, setup_branched_reaction_tree
):
    tree = setup_branched_reaction_tree()
    calc_patch = mocker.patch("aizynthfinder.context.scoring.scorers_mols.SCScore")
    calc_patch.return_value.return_value = 1.0

    scorer = DeltaSyntheticComplexityScorer(default_config, "dummy")

    assert pytest.approx(scorer(tree), abs=1e-3) == 0.2727


def test_delta_complexity_scorer_node(default_config, mocker, setup_branched_mcts):
    _, node = setup_branched_mcts()
    calc_patch = mocker.patch("aizynthfinder.context.scoring.scorers_mols.SCScore")
    calc_patch.return_value.return_value = 1.0

    scorer = DeltaSyntheticComplexityScorer(default_config, "dummy")

    assert pytest.approx(scorer(node), abs=1e-3) == 0.2727


def test_reaction_class_rank(default_config, load_reaction_tree, tmpdir):
    # Only test with a tree for this scorer as it is the same
    # functionality for node. Also don't need extensive testing,
    # because it is a wrapper around rxnutils.
    class_rank_path = str(tmpdir / "ranks.csv")
    with open(class_rank_path, "w") as fileobj:
        fileobj.write("reaction_class,rank_score\n")
        fileobj.write("abc,6\n")
        fileobj.write("xyz,5\n")

    preferred_classes_path = str(tmpdir / "classes.txt")
    with open(preferred_classes_path, "w") as fileobj:
        fileobj.write("abc\n")

    tree = ReactionTree.from_dict(
        load_reaction_tree("linear_route_w_metadata.json", remove_metadata=False)
    )

    scorer = ReactionClassRankScorer(
        default_config, class_rank_path, preferred_classes_path
    )

    assert scorer(tree) == pytest.approx(0.6666, abs=1e-4)


def test_deepset_scorer(default_config, load_reaction_tree, tmpdir, mocker):
    # Only test with a tree for this scorer as it is the same
    # functionality for node. Also don't need extensive testing,
    # because it is a wrapper around rxnutils.
    class_rank_path = str(tmpdir / "ranks.csv")
    with open(class_rank_path, "w") as fileobj:
        fileobj.write("reaction_class,rank_score\n")
        fileobj.write("0.0,6\n")

    tree = ReactionTree.from_dict(
        load_reaction_tree("linear_route_w_metadata.json", remove_metadata=False)
    )
    # The classifications in the json is unfortunately note valid classifications
    for rxn in tree.reactions():
        rxn.metadata["classification"] = "0.0"

    mocked_onxx_model = mocker.patch(
        "rxnutils.routes.deepset.scoring.onnxruntime.InferenceSession"
    )
    # The first return is for the SCScore, the second for the DeepSet model
    mocked_onxx_model.return_value.run.side_effect = [[[2]], [5.0]]
    scorer = DeepSetScorer(
        default_config,
        "dummy",
        "dummy",
        class_rank_path,
    )

    assert scorer(tree) == pytest.approx(3.99, abs=1e-2)
