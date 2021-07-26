import os
from tarfile import TarFile

import numpy as np
import pytest

from aizynthfinder.analysis import TreeAnalysis, RouteCollection
from aizynthfinder.analysis.utils import RouteSelectionArguments
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsSearchTree
from aizynthfinder.context.scoring import StateScorer, NumberOfReactionsScorer


def test_sort_nodes(setup_analysis):
    analysis, nodes = setup_analysis()

    best_nodes, best_scores = analysis.sort()

    assert len(best_nodes) == 7
    assert np.round(best_scores[0], 4) == 0.9940

    best_nodes, best_scores = analysis.sort(RouteSelectionArguments(nmin=0))

    assert len(best_nodes) == 0


def test_sort_nodes_nmax(setup_analysis):
    analysis, nodes = setup_analysis()

    best_nodes, best_scores = analysis.sort(RouteSelectionArguments(nmax=5))

    assert len(best_nodes) == 5


def test_sort_nodes_return_all(setup_analysis):
    analysis, nodes = setup_analysis()

    best_nodes, best_scores = analysis.sort(RouteSelectionArguments(return_all=True))

    assert len(best_nodes) == 1
    assert np.round(best_scores[0], 4) == 0.9940


def test_sort_nodes_scorer(setup_analysis):
    analysis, _ = setup_analysis(scorer=NumberOfReactionsScorer())

    best_nodes, best_scores = analysis.sort()

    assert len(best_nodes) == 10
    assert best_scores[0] == 2


def test_sort_routes(setup_analysis_andor_tree):
    analysis = setup_analysis_andor_tree(scorer=NumberOfReactionsScorer())

    best_routes, best_scores = analysis.sort()

    assert len(best_routes) == 3
    assert best_scores[0] == 1


def test_best_node(setup_analysis):
    analysis, nodes = setup_analysis(scorer=NumberOfReactionsScorer())

    assert analysis.best() is nodes[51]


def test_best_route(setup_analysis_andor_tree):
    analysis = setup_analysis_andor_tree(scorer=NumberOfReactionsScorer())

    assert len(list(analysis.best().reactions())) == 1


def test_tree_statistics(setup_analysis):
    analysis, _ = setup_analysis(scorer=NumberOfReactionsScorer())

    stats = analysis.tree_statistics()

    assert stats["number_of_nodes"] == 61
    assert stats["max_transforms"] == 7
    assert stats["max_children"] == 10
    assert stats["number_of_routes"] == 10
    assert stats["number_of_solved_routes"] == 1
    assert stats["top_score"] == 2
    assert stats["is_solved"]
    assert stats["number_of_steps"] == 2
    assert stats["number_of_precursors"] == 4
    assert stats["number_of_precursors_in_stock"] == 4
    mol_str = "CN1CCC(Cl)CC1, O, N#Cc1cccc(N)c1F, O=C(Cl)c1ccc(F)cc1"
    assert stats["precursors_in_stock"] == mol_str
    assert stats["precursors_not_in_stock"] == ""
    assert stats["policy_used_counts"] == {}


def test_tree_statistics_andor_tree(setup_analysis_andor_tree):
    analysis = setup_analysis_andor_tree(scorer=NumberOfReactionsScorer())

    stats = analysis.tree_statistics()
    assert stats["number_of_nodes"] == 9
    assert stats["max_transforms"] == 2
    assert stats["max_children"] == 2
    assert stats["number_of_routes"] == 3
    assert stats["number_of_solved_routes"] == 3
    assert stats["top_score"] == 1
    assert stats["is_solved"]
    assert stats["number_of_steps"] == 1
    assert stats["number_of_precursors"] == 2
    assert stats["number_of_precursors_in_stock"] == 2
    mol_str = "Cc1ccc2nc3ccccc3c(Cl)c2c1, Nc1ccc(NC(=S)Nc2ccccc2)cc1"
    assert stats["precursors_in_stock"] == mol_str
    assert stats["precursors_not_in_stock"] == ""
    assert stats["policy_used_counts"] == {}


def test_create_route_collection_full(setup_analysis, mocker):
    analysis, _ = setup_analysis()

    routes = RouteCollection.from_analysis(analysis)

    assert len(routes) == 7
    # Check a few of the routes
    assert np.round(routes.scores[0], 3) == 0.994
    assert len(routes.reaction_trees[0].graph) == 8
    assert np.round(routes.scores[1], 3) == 0.681
    assert len(routes.reaction_trees[1].graph) == 5

    assert "dict" not in routes[0]
    assert "json" not in routes[0]
    assert "image" not in routes[0]

    mocker.patch("aizynthfinder.reactiontree.ReactionTree.to_dict")
    mocker.patch("aizynthfinder.reactiontree.json.dumps")
    mocker.patch("aizynthfinder.utils.image.make_graphviz_image")

    # Just see that the code does not crash, does not verify content
    assert len(routes.images) == 7
    assert len(routes.dicts) == 7
    assert len(routes.jsons) == 7


def test_create_route_collection_andor_tree(setup_analysis_andor_tree):
    analysis = setup_analysis_andor_tree()

    routes = RouteCollection.from_analysis(analysis)

    assert len(routes) == 3
    assert routes.nodes == [None, None, None]


def test_compute_new_score(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis)

    routes.compute_scores(NumberOfReactionsScorer())

    assert np.round(routes.scores[0], 3) == 0.994
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.994
    assert routes.all_scores[0]["number of reactions"] == 2

    assert np.round(routes.scores[1], 3) == 0.681
    assert np.round(routes.all_scores[1]["state score"], 3) == 0.681
    assert routes.all_scores[1]["number of reactions"] == 1


def test_rescore_collection(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis)

    routes.rescore(NumberOfReactionsScorer())

    assert routes.scores[0] == 1
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.681
    assert routes.all_scores[0]["number of reactions"] == 1

    assert np.round(routes.all_scores[1]["state score"], 3) == 0.523
    assert routes.scores[1] == 1
    assert routes.all_scores[1]["number of reactions"] == 1


def test_dict_with_scores(setup_analysis):
    analysis, _ = setup_analysis()
    routes = RouteCollection.from_analysis(analysis)

    dicts = routes.dict_with_scores()

    assert "scores" not in routes.dicts[0]
    assert "scores" in dicts[0]
    assert np.round(dicts[0]["scores"]["state score"], 3) == 0.994


def test_compute_new_score_for_trees(default_config, setup_linear_reaction_tree):
    rt = setup_linear_reaction_tree()
    routes = RouteCollection(reaction_trees=[rt])

    assert routes.nodes[0] is None
    assert routes.scores[0] is np.nan
    assert routes.all_scores[0] == {}

    routes.compute_scores(StateScorer(default_config), NumberOfReactionsScorer())

    assert routes.scores[0] is np.nan
    assert np.round(routes.all_scores[0]["state score"], 3) == 0.994
    assert routes.all_scores[0]["number of reactions"] == 2


def test_rescore_collection_for_trees(default_config, setup_linear_reaction_tree):
    rt = setup_linear_reaction_tree()
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
    assert len(combined_dict["children"][1]["children"][1]["children"]) == 2
    assert combined_dict["children"][1]["children"][1]["children"][0]["is_reaction"]
    assert combined_dict == expected


def test_create_combine_tree_dict_from_tree(
    setup_stock, default_config, load_reaction_tree, shared_datadir
):
    setup_stock(
        default_config,
        "Nc1ccc(NC(=S)Nc2ccccc2)cc1",
        "Cc1ccc2nc3ccccc3c(Cl)c2c1",
        "Nc1ccc(N)cc1",
        "S=C=Nc1ccccc1",
        "Cc1ccc2nc3ccccc3c(N)c2c1",
        "Nc1ccc(Br)cc1",
    )
    search_tree = MctsSearchTree.from_json(
        shared_datadir / "tree_for_clustering.json", default_config
    )
    analysis = TreeAnalysis(search_tree)
    collection = RouteCollection.from_analysis(
        analysis, RouteSelectionArguments(nmin=3)
    )
    expected = load_reaction_tree("combined_example_tree.json")

    combined_dict = collection.combined_reaction_trees().to_dict()

    assert len(combined_dict["children"]) == 2
    assert combined_dict["children"][0]["is_reaction"]
    assert len(combined_dict["children"][0]["children"]) == 2
    assert len(combined_dict["children"][1]["children"]) == 2
    assert len(combined_dict["children"][1]["children"][1]["children"]) == 2
    assert combined_dict["children"][1]["children"][1]["children"][0]["is_reaction"]
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


def test_clustering_collection_timeout(load_reaction_tree):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(
                load_reaction_tree("routes_for_clustering.json", idx)
            )
            for idx in range(3)
        ]
    )
    cluster_labels = collection.cluster(n_clusters=1, timeout=0)

    assert len(cluster_labels) == 0
    assert collection.clusters is None
