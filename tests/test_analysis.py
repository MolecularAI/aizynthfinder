import json

import numpy as np

from aizynthfinder.chem import Molecule
from aizynthfinder.analysis import TreeAnalysis, ReactionTree, RouteCollection
from aizynthfinder.mcts.mcts import SearchTree
from aizynthfinder.scoring import StateScorer, NumberOfReactionsScorer


def test_sort_nodes(setup_analysis):
    analysis, nodes = setup_analysis()

    best_nodes, best_scores = analysis.sort_nodes()

    assert len(best_nodes) == 7
    assert np.round(best_scores[0], 4) == 0.9940

    best_nodes, best_scores = analysis.sort_nodes(min_return=0)

    assert len(best_nodes) == 0


def test_sort_nodes_scorer(setup_analysis):
    analysis, _ = setup_analysis(scorer=NumberOfReactionsScorer())

    best_nodes, best_scores = analysis.sort_nodes()

    assert len(best_nodes) == 10
    assert best_scores[0] == 2


def test_best_node(setup_analysis):
    analysis, nodes = setup_analysis(scorer=NumberOfReactionsScorer())

    assert analysis.best_node() is nodes[51]


def test_tree_statistics(setup_analysis):
    analysis, _ = setup_analysis(scorer=NumberOfReactionsScorer())

    stats = analysis.tree_statistics()

    assert stats["number_of_nodes"] == 61
    assert stats["max_transforms"] == 7
    assert stats["max_children"] == 10
    assert stats["number_of_leafs"] == 10
    assert stats["number_of_solved_leafs"] == 1
    assert stats["top_score"] == 2
    assert stats["is_solved"]
    assert stats["number_of_steps"] == 2
    assert stats["number_of_precursors"] == 4
    assert stats["number_of_precursors_in_stock"] == 4
    mol_str = "CN1CCC(Cl)CC1, O, N#Cc1cccc(N)c1F, O=C(Cl)c1ccc(F)cc1"
    assert stats["precursors_in_stock"] == mol_str
    assert stats["precursors_not_in_stock"] == ""
    assert stats["policy_used_counts"] == {}


def test_route_to_reactiontree(setup_analysis):
    analysis, _ = setup_analysis()

    reaction_tree = ReactionTree.from_analysis(analysis)

    mol_nodes = list(reaction_tree.molecules())
    assert len(mol_nodes) == 6

    reaction_nodes = list(reaction_tree.reactions())
    assert len(reaction_nodes) == 2

    leaf_nodes = list(reaction_tree.leafs())
    assert len(leaf_nodes) == 4
    assert set(node.smiles for node in leaf_nodes) == {
        "N#Cc1cccc(N)c1F",
        "O=C(Cl)c1ccc(F)cc1",
        "CN1CCC(Cl)CC1",
        "O",
    }


def test_reactiontree_to_json(setup_analysis, load_reaction_tree):
    expected = load_reaction_tree("sample_reaction.json")
    analysis, _ = setup_analysis()

    resp = ReactionTree.from_analysis(analysis).to_json()
    assert json.loads(resp) == expected


def test_reactiontree_from_dict(load_reaction_tree):
    expected = load_reaction_tree("sample_reaction.json")

    rt = ReactionTree.from_dict(expected)

    # Simply check that the to_dict() and from_dict() gives/produces the same dict
    resp = rt.to_dict()
    assert resp == expected


def test_reactiontree_to_image(load_reaction_tree, mocker):
    patched_add_mol = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_molecule"
    )
    patched_add_reaction = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_reaction"
    )
    patched_add_edge = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_edge"
    )
    mocker.patch("aizynthfinder.utils.image.GraphvizReactionGraph.to_image")

    tree = load_reaction_tree("sample_reaction.json")
    rt = ReactionTree.from_dict(tree)

    rt.to_image()

    assert patched_add_mol.call_count == len(list(rt.molecules()))
    assert patched_add_reaction.call_count == len(list(rt.reactions()))
    assert patched_add_edge.call_count == len(rt.graph.edges)


def test_reactiontree_to_image_hiding(load_reaction_tree, mocker):
    patched_add_mol = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_molecule"
    )
    patched_add_reaction = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_reaction"
    )
    patched_add_edge = mocker.patch(
        "aizynthfinder.utils.image.GraphvizReactionGraph.add_edge"
    )
    mocker.patch("aizynthfinder.utils.image.GraphvizReactionGraph.to_image")

    tree = load_reaction_tree("finder_output_mol2.json", 1)
    rt = ReactionTree.from_dict(tree)
    assert rt.has_repeating_patterns

    rt.to_image(show_all=True)

    assert patched_add_mol.call_count == len(list(rt.molecules()))
    assert patched_add_reaction.call_count == len(list(rt.reactions()))
    assert patched_add_edge.call_count == len(rt.graph.edges)

    patched_add_mol.reset_mock()
    patched_add_reaction.reset_mock()
    patched_add_edge.reset_mock()

    rt.to_image(show_all=False)

    assert patched_add_mol.call_count == len(list(rt.molecules())) - 3
    assert patched_add_reaction.call_count == len(list(rt.reactions())) - 2
    assert patched_add_edge.call_count == len(rt.graph.edges) - 5


def test_find_repetetive_patterns(load_reaction_tree):
    tree_with_repetetive_patterns = load_reaction_tree("finder_output_mol2.json", 1)

    rt = ReactionTree.from_dict(tree_with_repetetive_patterns)

    assert rt.has_repeating_patterns
    assert len([node for node in rt.graph if rt.graph.nodes[node]["hide"]]) == 5


def test_find_repetetive_patterns_no_patterns(load_reaction_tree):
    tree_with_no_repetetive_patterns = load_reaction_tree("finder_output_mol2.json", 0)

    rt = ReactionTree.from_dict(tree_with_no_repetetive_patterns)

    assert not rt.has_repeating_patterns
    assert len([node for node in rt.graph if rt.graph.nodes[node]["hide"]]) == 0


def test_find_repetetive_patterns_created_tree(
    default_config, mock_stock, shared_datadir
):
    mock_stock(default_config, Molecule(smiles="CC"), Molecule(smiles="C"))

    # Try one with 2 repetetive units
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_with_repetition.json", default_config
    )
    analysis = TreeAnalysis(search_tree)

    rt = ReactionTree.from_analysis(analysis)

    assert rt.has_repeating_patterns
    hidden_nodes = [
        node for node in rt.graph if rt.graph.nodes[node].get("hide", False)
    ]
    assert len(hidden_nodes) == 5

    # Try one with 3 repetetive units
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_with_3_repetitions.json", default_config
    )
    analysis = TreeAnalysis(search_tree)

    rt = ReactionTree.from_analysis(analysis)

    assert rt.has_repeating_patterns
    hidden_nodes = [
        node for node in rt.graph if rt.graph.nodes[node].get("hide", False)
    ]
    assert len(hidden_nodes) == 10


def test_find_repetetive_patterns_created_tree_no_patterns(
    default_config, mock_stock, shared_datadir
):
    mock_stock(default_config, Molecule(smiles="CC"), Molecule(smiles="CCCO"))

    # Try with a short tree (3 nodes, 1 reaction)
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    analysis = TreeAnalysis(search_tree)

    rt = ReactionTree.from_analysis(analysis)

    assert not rt.has_repeating_patterns
    hidden_nodes = [
        node for node in rt.graph if rt.graph.nodes[node].get("hide", False)
    ]
    assert len(hidden_nodes) == 0

    # Try with something longer
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition_longer.json", default_config
    )
    analysis = TreeAnalysis(search_tree)

    rt = ReactionTree.from_analysis(analysis)

    assert not rt.has_repeating_patterns


def test_route_node_depth(load_reaction_tree):
    dict_ = load_reaction_tree("finder_output_mol2.json", 0)
    rt = ReactionTree.from_dict(dict_)

    mols = list(rt.molecules())

    assert rt.depth(mols[0]) == 0
    assert rt.depth(mols[1]) == 2
    assert rt.depth(mols[2]) == 4
    assert rt.depth(mols[3]) == 6
    assert rt.depth(mols[4]) == 6
    assert rt.depth(mols[5]) == 8
    assert rt.depth(mols[6]) == 8
    assert rt.depth(mols[7]) == 10
    assert rt.depth(mols[8]) == 12
    assert rt.depth(mols[9]) == 12
    assert rt.depth(mols[10]) == 2


def test_create_route_collection_full(setup_analysis, mocker):
    analysis, _ = setup_analysis()

    routes = RouteCollection.from_analysis(analysis, 5)

    assert len(routes) == 7
    # Check a few of the routes
    assert np.round(routes.scores[0], 3) == 0.994
    assert len(routes.reaction_trees[0].graph) == 8
    assert np.round(routes.scores[1], 3) == 0.681
    assert len(routes.reaction_trees[1].graph) == 5

    assert "dict" not in routes[0]
    assert "json" not in routes[0]
    assert "image" not in routes[0]

    mocker.patch("aizynthfinder.analysis.ReactionTree.to_dict")
    mocker.patch("aizynthfinder.analysis.json.dumps")
    mocker.patch("aizynthfinder.utils.image.GraphvizReactionGraph.to_image")
    mocker.patch("aizynthfinder.utils.image.GraphvizReactionGraph.add_molecule")

    # Just see that the code does not crash, does not verify content
    assert len(routes.images) == 7
    assert len(routes.dicts) == 7
    assert len(routes.jsons) == 7


def test_compute_new_score(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis, 5)

    routes.compute_scores(NumberOfReactionsScorer())

    assert np.round(routes.scores[0], 3) == 0.994
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.994
    assert routes.all_scores[0]["number of reactions"] == 2

    assert np.round(routes.scores[1], 3) == 0.681
    assert np.round(routes.all_scores[1]["state score"], 3) == 0.681
    assert routes.all_scores[1]["number of reactions"] == 1


def test_rescore_collection(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis, 5)

    routes.rescore(NumberOfReactionsScorer())

    assert routes.scores[0] == 1
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.681
    assert routes.all_scores[0]["number of reactions"] == 1

    assert np.round(routes.all_scores[1]["state score"], 3) == 0.523
    assert routes.scores[1] == 1
    assert routes.all_scores[1]["number of reactions"] == 1


def test_dict_with_scores(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis, 5)

    dicts = routes.dict_with_scores()

    assert "scores" not in routes.dicts[0]
    assert "scores" in dicts[0]
    assert np.round(dicts[0]["scores"]["state score"], 3) == 0.994


def test_compute_new_score_for_trees(default_config, mock_stock, load_reaction_tree):
    mock_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )
    rt = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    routes = RouteCollection(reaction_trees=[rt])

    assert routes.nodes[0] is None
    assert routes.scores[0] is np.nan
    assert routes.all_scores[0] == {}

    routes.compute_scores(StateScorer(default_config), NumberOfReactionsScorer())

    assert routes.scores[0] is np.nan
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.994
    assert routes.all_scores[0]["number of reactions"] == 2


def test_rescore_collection_for_trees(default_config, mock_stock, load_reaction_tree):
    mock_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )
    rt = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    routes = RouteCollection(reaction_trees=[rt])
    routes.compute_scores(StateScorer(default_config))

    routes.rescore(NumberOfReactionsScorer())

    assert routes.scores[0] == 2
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.994
    assert routes.all_scores[0]["number of reactions"] == 2
