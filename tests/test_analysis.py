import json
import os
from tarfile import TarFile

import numpy as np
import pytest

from aizynthfinder.chem import Molecule
from aizynthfinder.analysis import TreeAnalysis, ReactionTree, RouteCollection
from aizynthfinder.mcts.mcts import SearchTree
from aizynthfinder.context.scoring import StateScorer, NumberOfReactionsScorer


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
    patched_make_image = mocker.patch("aizynthfinder.analysis.make_graphviz_image")

    tree = load_reaction_tree("sample_reaction.json")
    rt = ReactionTree.from_dict(tree)

    rt.to_image()

    patched_make_image.assert_called_once()
    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules()))
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions()))


def test_reactiontree_to_image_hiding(load_reaction_tree, mocker):
    patched_make_image = mocker.patch("aizynthfinder.analysis.make_graphviz_image")

    tree = load_reaction_tree("sample_reaction_with_hidden.json", 1)
    rt = ReactionTree.from_dict(tree)
    assert rt.has_repeating_patterns

    rt.to_image(show_all=True)

    patched_make_image.assert_called_once()
    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules()))
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions()))

    patched_make_image.reset_mock()

    rt.to_image(show_all=False)

    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules())) - 3
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions())) - 2


def test_find_repetetive_patterns(load_reaction_tree):
    tree_with_repetetive_patterns = load_reaction_tree(
        "sample_reaction_with_hidden.json", 1
    )

    rt = ReactionTree.from_dict(tree_with_repetetive_patterns)

    assert rt.has_repeating_patterns
    assert len([node for node in rt.graph if rt.graph.nodes[node]["hide"]]) == 5


def test_find_repetetive_patterns_no_patterns(load_reaction_tree):
    tree_with_no_repetetive_patterns = load_reaction_tree(
        "sample_reaction_with_hidden.json", 0
    )

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


def test_route_node_depth_from_analysis(default_config, mock_stock, shared_datadir):
    mock_stock(default_config, Molecule(smiles="CC"), Molecule(smiles="CCCO"))
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    analysis = TreeAnalysis(search_tree)
    rt = ReactionTree.from_analysis(analysis)

    mols = list(rt.molecules())

    assert rt.depth(mols[0]) == 0
    assert rt.depth(mols[1]) == 2
    assert rt.depth(mols[2]) == 2

    rxns = list(rt.reactions())
    assert rt.depth(rxns[0]) == 1

    for mol in rt.molecules():
        assert rt.depth(mol) == 2 * rt.graph.nodes[mol]["transform"]


def test_route_node_depth_from_json(load_reaction_tree):
    dict_ = load_reaction_tree("sample_reaction_with_hidden.json", 0)
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

    rxns = list(rt.reactions())

    assert rt.depth(rxns[0]) == 1
    assert rt.depth(rxns[1]) == 3

    for mol in rt.molecules():
        assert rt.depth(mol) == 2 * rt.graph.nodes[mol]["transform"]


def test_route_distance_self(load_reaction_tree):
    dict_ = load_reaction_tree("sample_reaction_with_hidden.json", 0)
    rt = ReactionTree.from_dict(dict_)

    assert rt.distance_to(rt) == 0.0


def test_route_distance_other(load_reaction_tree):
    dict_ = load_reaction_tree("routes_for_clustering.json", 0)
    rt1 = ReactionTree.from_dict(dict_)
    dict_ = load_reaction_tree("routes_for_clustering.json", 1)
    rt2 = ReactionTree.from_dict(dict_)

    dist = rt1.distance_to(rt2, content="molecules")

    assert pytest.approx(dist, abs=1e-2) == 2.6522


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
    mocker.patch("aizynthfinder.utils.image.make_graphviz_image")

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


def test_create_combine_tree_dict_from_json(load_reaction_tree):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0)),
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 1)),
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 2)),
        ]
    )
    expected = load_reaction_tree("combined_example_tree.json")

    combined_dict = collection.combined_reaction_trees().to_dict()

    assert len(combined_dict["children"]) == 2
    assert combined_dict["children"][0]["is_reaction"]
    assert len(combined_dict["children"][0]["children"]) == 2
    assert len(combined_dict["children"][1]["children"]) == 2
    assert len(combined_dict["children"][1]["children"][0]["children"]) == 2
    assert combined_dict["children"][1]["children"][0]["children"][0]["is_reaction"]
    assert combined_dict == expected


def test_create_combine_tree_dict_from_tree(
    mock_stock, default_config, load_reaction_tree, shared_datadir
):
    mock_stock(
        default_config,
        "Nc1ccc(NC(=S)Nc2ccccc2)cc1",
        "Cc1ccc2nc3ccccc3c(Cl)c2c1",
        "Nc1ccc(N)cc1",
        "S=C=Nc1ccccc1",
        "Cc1ccc2nc3ccccc3c(N)c2c1",
        "Nc1ccc(Br)cc1",
    )
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_for_clustering.json", default_config
    )
    analysis = TreeAnalysis(search_tree)
    collection = RouteCollection.from_analysis(analysis, 3)
    expected = load_reaction_tree("combined_example_tree.json")

    combined_dict = collection.combined_reaction_trees().to_dict()

    assert len(combined_dict["children"]) == 2
    assert combined_dict["children"][0]["is_reaction"]
    assert len(combined_dict["children"][0]["children"]) == 2
    assert len(combined_dict["children"][1]["children"]) == 2
    assert len(combined_dict["children"][1]["children"][0]["children"]) == 2
    assert combined_dict["children"][1]["children"][0]["children"][0]["is_reaction"]
    assert combined_dict == expected


def test_create_combine_tree_to_visjs(load_reaction_tree, tmpdir):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0)),
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 1)),
            ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 2)),
        ]
    )
    tar_filename = str(tmpdir / "routes.tar")
    combined = collection.combined_reaction_trees()

    combined.to_visjs_page(tar_filename)

    assert os.path.exists(tar_filename)
    with TarFile(tar_filename) as tarobj:
        assert "./route.html" in tarobj.getnames()
        assert len([name for name in tarobj.getnames() if name.endswith(".png")]) == 8


def test_distance_collection(load_reaction_tree):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(
                load_reaction_tree("routes_for_clustering.json", idx)
            )
            for idx in range(3)
        ]
    )

    dist_mat1 = collection.distance_matrix()
    dist_mat2 = collection.distance_matrix(recreate=True)

    assert (dist_mat1 - dist_mat2).sum() == 0

    dist_mat3 = collection.distance_matrix(content="molecules")

    assert (dist_mat1 - dist_mat3).sum() != 0
    assert len(dist_mat3) == 3
    assert pytest.approx(dist_mat3[0, 1], abs=1e-2) == 2.6522
    assert pytest.approx(dist_mat3[0, 2], abs=1e-2) == 3.0779
    assert pytest.approx(dist_mat3[2, 1], abs=1e-2) == 0.7483


def test_clustering_collection(load_reaction_tree):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(
                load_reaction_tree("routes_for_clustering.json", idx)
            )
            for idx in range(3)
        ]
    )
    collection.cluster(n_clusters=1)

    assert len(collection.clusters) == 2
    assert collection.clusters[0].reaction_trees == collection.reaction_trees[1:3]
    assert collection.clusters[1].reaction_trees == [collection.reaction_trees[0]]
