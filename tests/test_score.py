import pytest
import numpy as np

from aizynthfinder.scoring import (
    StateScorer,
    NumberOfReactionsScorer,
    AverageTemplateOccurenceScorer,
    ScorerCollection,
    ScorerException,
)
from aizynthfinder.chem import Molecule
from aizynthfinder.mcts.mcts import SearchTree
from aizynthfinder.analysis import ReactionTree


def test_state_scorer_node(generate_root):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    scorer = StateScorer()

    assert repr(scorer) == "state score"
    assert round(scorer(root), 4) == 0.0491


def test_state_scorer_tree(load_reaction_tree, default_config, mock_stock):
    mock_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = StateScorer(default_config)

    assert round(scorer(tree), 4) == 0.994


def test_sort(shared_datadir, default_config, mock_stock):
    mock_stock(default_config, "CCCO", "CC")
    search_tree = SearchTree.from_json(
        shared_datadir / "tree_without_repetition.json", default_config
    )
    nodes = list(search_tree.graph())
    scorer = StateScorer()

    sorted_nodes, scores = scorer.sort(nodes)

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

    assert scorer(nodes[1]) == 5


def test_template_occurence_scorer_tree(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("sample_reaction.json"))
    scorer = AverageTemplateOccurenceScorer()

    assert scorer(tree) == 0


def test_create_scorer_collection(default_config):
    collection = ScorerCollection(default_config)

    assert len(collection) == 3

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

    collection.add(StateScorer())

    assert "state score" in collection.names()


def test_add_scorer_to_collection_no_scorer(default_config):
    collection = ScorerCollection(default_config)

    with pytest.raises(ScorerException):
        collection.add(Molecule(smiles="CCC"))


def test_load_scorer_to_collection_only_class(default_config):
    collection = ScorerCollection(default_config)
    del collection["state score"]

    collection.load(**{"StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_full_package(default_config):
    collection = ScorerCollection(default_config)
    del collection["state score"]

    collection.load(**{"aizynthfinder.scoring.StateScorer": {}})

    assert "state score" in collection.names()


def test_load_scorer_to_collection_failures(default_config):
    collection = ScorerCollection(default_config)

    with pytest.raises(ScorerException, match=".*load module.*"):
        collection.load(**{"mypackage.scoring.StateScorer": {}})

    with pytest.raises(ScorerException, match=".*class.*"):
        collection.load(**{"aizynthfinder.scoring.NoScorer": {}})
