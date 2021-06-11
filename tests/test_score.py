import pytest
import numpy as np

from aizynthfinder.context.scoring import (
    StateScorer,
    NumberOfReactionsScorer,
    AverageTemplateOccurenceScorer,
    NumberOfPrecursorsScorer,
    NumberOfPrecursorsInStockScorer,
    PriceSumScorer,
    RouteCostScorer,
    ScorerCollection,
    ScorerException,
)
from aizynthfinder.chem import Molecule, UniqueMolecule
from aizynthfinder.mcts.mcts import SearchTree
from aizynthfinder.reactiontree import ReactionTree


def test_state_scorer_node(generate_root, default_config):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    scorer = StateScorer(default_config)

    assert repr(scorer) == "state score"
    assert round(scorer(root), 4) == 0.0491


def test_state_scorer_nodes(generate_root, default_config):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    scorer = StateScorer(default_config)

    scores = scorer([root, root])

    assert repr(scorer) == "state score"
    assert round(scores[0], 4) == 0.0491
    assert round(scores[1], 4) == 0.0491


def test_state_scorer_tree(load_reaction_tree, default_config, mock_stock):
    mock_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = StateScorer(default_config)

    assert round(scorer(tree), 4) == 0.994


def test_state_scorer_trees(load_reaction_tree, default_config, mock_stock):
    mock_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = StateScorer(default_config)

    scores = scorer([tree, tree])

    assert round(scores[0], 4) == 0.994
    assert round(scores[1], 4) == 0.994


def test_sort(shared_datadir, default_config, mock_stock):
    mock_stock(default_config, "CCCO", "CC")
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    nodes = list(search_tree.graph())
    scorer = StateScorer(default_config)

    sorted_nodes, scores, _ = scorer.sort(nodes)

    assert [np.round(score, 4) for score in scores] == [0.9976, 0.0491]
    assert sorted_nodes == [nodes[1], nodes[0]]


def test_number_of_reaction_scorer_node(shared_datadir, default_config):
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    nodes = list(search_tree.graph())
    scorer = NumberOfReactionsScorer()

    assert scorer(nodes[1]) == 1


def test_number_of_reaction_scorer_tree(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = NumberOfReactionsScorer()

    assert scorer(tree) == 2


def test_template_occurence_scorer_no_metadata(shared_datadir, default_config):
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    nodes = list(search_tree.graph())
    scorer = AverageTemplateOccurenceScorer()

    assert scorer(nodes[1]) == 0


def test_template_occurence_scorer(shared_datadir, default_config):
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    nodes = list(search_tree.graph())
    nodes[0][nodes[1]]["action"].metadata["library_occurence"] = 5
    scorer = AverageTemplateOccurenceScorer()

    assert scorer(nodes[0]) == 0
    assert scorer(nodes[1]) == 5


def test_template_occurence_scorer_tree(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = AverageTemplateOccurenceScorer()

    assert scorer(tree) == 0


def test_template_occurence_scorer_tree_one_node():
    rt = ReactionTree()
    rt.root = Molecule(smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    rt.graph.add_node(rt.root)
    scorer = AverageTemplateOccurenceScorer()

    assert scorer(rt) == 0.0


def test_scorers_one_mcts_node(default_config):
    tree = SearchTree(default_config, root_smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    node = tree.root

    assert pytest.approx(StateScorer(default_config)(node), abs=1e-3) == 0.0497
    assert NumberOfReactionsScorer(default_config)(node) == 0
    assert NumberOfPrecursorsScorer(default_config)(node) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(node) == 0
    assert PriceSumScorer(default_config)(node) == 10
    assert RouteCostScorer(default_config)(node) == 10


def test_scoring_branched_mcts_tree(shared_datadir, default_config):
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_with_branching.json", default_config
    )
    nodes = list(search_tree.graph())

    assert pytest.approx(StateScorer(default_config)(nodes[-1]), abs=1e-6) == 0.00012363
    assert NumberOfReactionsScorer()(nodes[-1]) == 14
    assert NumberOfPrecursorsScorer(default_config)(nodes[-1]) == 8
    assert NumberOfPrecursorsInStockScorer(default_config)(nodes[-1]) == 0
    assert PriceSumScorer(default_config)(nodes[-1]) == 80
    cost_score = RouteCostScorer(default_config)(nodes[-1])
    assert pytest.approx(cost_score, abs=1e-3) == 410.6577


def test_scoring_branch_mcts_tree_in_stock(shared_datadir, default_config, mock_stock):
    mock_stock(
        default_config,
        "CC(C)(C)CO",
        "CC(C)(C)OC(=O)N(CCCl)CCCl",
        "N#CCc1cccc(O)c1F",
        "O=[N+]([O-])c1ccccc1F",
        "O=C1CCC(=O)N1Br",
        "O=C=Nc1csc(C(F)(F)F)n1",
        "CCC[Sn](Cl)(CCC)CCC",
        "COc1ccc2ncsc2c1",
    )
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_with_branching.json", default_config
    )
    nodes = list(search_tree.graph())

    assert pytest.approx(StateScorer(default_config)(nodes[-1]), abs=1e-3) == 0.950
    assert NumberOfReactionsScorer()(nodes[-1]) == 14
    assert NumberOfPrecursorsScorer(default_config)(nodes[-1]) == 8
    assert NumberOfPrecursorsInStockScorer(default_config)(nodes[-1]) == 8
    assert PriceSumScorer(default_config)(nodes[-1]) == 8
    cost_score = RouteCostScorer(default_config)(nodes[-1])
    assert pytest.approx(cost_score, abs=1e-3) == 77.4797


def test_scorers_tree_one_node_route(default_config):
    tree = ReactionTree()
    tree.root = UniqueMolecule(smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    tree.graph.add_node(tree.root)

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-3) == 0.0497
    assert NumberOfReactionsScorer(default_config)(tree) == 0
    assert NumberOfPrecursorsScorer(default_config)(tree) == 1
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 0
    assert PriceSumScorer(default_config)(tree) == 10
    assert RouteCostScorer(default_config)(tree) == 10


def test_scoring_branched_route(load_reaction_tree, default_config):
    tree = ReactionTree.from_dict(load_reaction_tree("branched_route.json"))

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-6) == 0.00012363
    assert NumberOfReactionsScorer(default_config)(tree) == 14
    assert NumberOfPrecursorsScorer(default_config)(tree) == 8
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 0
    assert PriceSumScorer(default_config)(tree) == 80
    cost_score = RouteCostScorer(default_config)(tree)
    assert pytest.approx(cost_score, abs=1e-3) == 410.6577


def test_scoring_branched_route_in_stock(
    load_reaction_tree, default_config, mock_stock
):
    mock_stock(
        default_config,
        "CC(C)(C)CO",
        "CC(C)(C)OC(=O)N(CCCl)CCCl",
        "N#CCc1cccc(O)c1F",
        "O=[N+]([O-])c1ccccc1F",
        "O=C1CCC(=O)N1Br",
        "O=C=Nc1csc(C(F)(F)F)n1",
        "CCC[Sn](Cl)(CCC)CCC",
        "COc1ccc2ncsc2c1",
    )
    tree = ReactionTree.from_dict(load_reaction_tree("branched_route.json"))

    assert pytest.approx(StateScorer(default_config)(tree), abs=1e-3) == 0.950
    assert NumberOfReactionsScorer(default_config)(tree) == 14
    assert NumberOfPrecursorsScorer(default_config)(tree) == 8
    assert NumberOfPrecursorsInStockScorer(default_config)(tree) == 8
    assert PriceSumScorer(default_config)(tree) == 8
    cost_score = RouteCostScorer(default_config)(tree)
    assert pytest.approx(cost_score, abs=1e-3) == 77.4797


def test_create_scorer_collection(default_config):
    collection = ScorerCollection(default_config)

    assert len(collection) == 5

    assert "state score" in collection.names()
    assert "number of reactions" in collection.names()

    assert isinstance(collection["state score"], StateScorer)

    with pytest.raises(KeyError):
        collection["dummy"]


def test_delete_scorer_to_collection(default_config):
    collection = ScorerCollection(default_config)

    del collection["state score"]

    assert "state score" not in collection.names()


def test_add_scorer_to_collection(default_config):
    collection = ScorerCollection(default_config)
    del collection["state score"]

    collection.load(StateScorer(default_config))

    assert "state score" in collection.names()


def test_add_scorer_to_collection_no_scorer(default_config):
    collection = ScorerCollection(default_config)

    with pytest.raises(ScorerException):
        collection.load(Molecule(smiles="CCC"))


def test_load_scorer_to_collection_only_class(default_config):
    collection = ScorerCollection(default_config)
    del collection["state score"]

    collection.load_from_config(**{"StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_full_package(default_config):
    collection = ScorerCollection(default_config)
    del collection["state score"]

    collection.load_from_config(**{"aizynthfinder.context.scoring.StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_failures(default_config):
    collection = ScorerCollection(default_config)

    with pytest.raises(ScorerException, match=".*load module.*"):
        collection.load_from_config(**{"mypackage.scoring.StateScorer": {}})

    with pytest.raises(ScorerException, match=".*class.*"):
        collection.load_from_config(**{"aizynthfinder.context.scoring.NoScorer": {}})
