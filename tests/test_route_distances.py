import pytest

from aizynthfinder.utils.route_clustering.distances import (
    AptedConfig,
    ReactionTreeWrapper,
    TreeContent,
)
from aizynthfinder.analysis import ReactionTree
from aizynthfinder.chem import Molecule


node1 = {"type": "mol", "fingerprint": [0, 1, 0], "children": ["A", "B", "C"]}

node2 = {"type": "mol", "fingerprint": [1, 1, 0]}


def test_rename_cost_different_types():
    config = AptedConfig()

    cost = config.rename({"type": "type1"}, {"type": "type2"})

    assert cost == 1


def test_rename_cost_same_types():
    config = AptedConfig()

    cost = config.rename(node1, node2)

    assert cost == 0.5


def test_get_children_fixed():
    config = AptedConfig()

    assert config.children(node1) == ["A", "B", "C"]


def test_get_children_random():
    config = AptedConfig(randomize=True)

    children = config.children(node1)

    assert len(children) == 3
    for expected_child in ["A", "B", "C"]:
        assert expected_child in children


@pytest.mark.parametrize(
    "route_index", [1, 2],
)
def test_create_wrapper(load_reaction_tree, route_index):
    tree = ReactionTree.from_dict(
        load_reaction_tree("routes_for_clustering.json", route_index)
    )

    wrapper = ReactionTreeWrapper(tree)

    assert wrapper.info["content"] == TreeContent.MOLECULES
    assert wrapper.info["tree count"] == 4
    assert wrapper.info["root"] is tree.root
    assert len(wrapper.trees) == 4

    wrapper = ReactionTreeWrapper(tree, TreeContent.REACTIONS)

    assert wrapper.info["content"] == TreeContent.REACTIONS
    assert wrapper.info["tree count"] == 1
    assert wrapper.info["root"] is list(tree.graph[tree.root])[0]
    assert len(wrapper.trees) == 1

    wrapper = ReactionTreeWrapper(tree, TreeContent.BOTH)

    assert wrapper.info["content"] == TreeContent.BOTH
    assert wrapper.info["tree count"] == 4
    assert wrapper.info["root"] is tree.root
    assert len(wrapper.trees) == 4


def test_create_wrapper_no_reaction():
    tree = ReactionTree()
    mol = Molecule(smiles="CCC")
    tree.graph.add_node(mol)
    tree.root = mol

    wrapper = ReactionTreeWrapper(tree)
    assert wrapper.info["tree count"] == 1
    assert wrapper.info["root"] is mol
    assert len(wrapper.trees) == 1

    wrapper = ReactionTreeWrapper(tree, TreeContent.REACTIONS)
    assert wrapper.info["tree count"] == 0
    assert wrapper.info["root"] is None
    assert len(wrapper.trees) == 0

    wrapper = ReactionTreeWrapper(tree, TreeContent.BOTH)
    assert wrapper.info["tree count"] == 1
    assert wrapper.info["root"] is mol
    assert len(wrapper.trees) == 1


def test_create_one_tree_of_molecules(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))

    wrapper = ReactionTreeWrapper(tree, exhaustive_limit=1)

    assert wrapper.info["tree count"] == 2
    assert len(wrapper.trees) == 1

    mol_nodes = list(tree.molecules())
    assert wrapper.first_tree["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.first_tree["children"]) == 2

    child_smiles = [child["smiles"] for child in wrapper.first_tree["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles


def test_create_one_tree_of_reactions(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))

    wrapper = ReactionTreeWrapper(
        tree, content=TreeContent.REACTIONS, exhaustive_limit=1
    )

    assert wrapper.info["tree count"] == 1
    assert len(wrapper.trees) == 1

    rxn_nodes = list(tree.reactions())
    assert wrapper.first_tree["smiles"] == rxn_nodes[0].smiles
    assert len(wrapper.first_tree["children"]) == 0


def test_create_one_tree_of_everything(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))

    wrapper = ReactionTreeWrapper(tree, content=TreeContent.BOTH, exhaustive_limit=1)

    assert wrapper.info["tree count"] == 2
    assert len(wrapper.trees) == 1

    mol_nodes = list(tree.molecules())
    rxn_nodes = list(tree.reactions())
    assert wrapper.first_tree["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.first_tree["children"]) == 1

    child1 = wrapper.first_tree["children"][0]
    assert child1["smiles"] == rxn_nodes[0].smiles
    assert len(child1["children"]) == 2

    child_smiles = [child["smiles"] for child in child1["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles


def test_create_all_trees_of_molecules(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))

    wrapper = ReactionTreeWrapper(tree)

    assert wrapper.info["tree count"] == 2
    assert len(wrapper.trees) == 2

    mol_nodes = list(tree.molecules())
    # Assert first tree
    assert wrapper.first_tree["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.first_tree["children"]) == 2

    child_smiles = [child["smiles"] for child in wrapper.first_tree["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles

    # Assert second tree
    assert wrapper.trees[1]["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.trees[1]["children"]) == 2

    child_smiles = [child["smiles"] for child in wrapper.trees[1]["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles[::-1]


def test_create_two_trees_of_everything(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))

    wrapper = ReactionTreeWrapper(tree, content=TreeContent.BOTH)

    assert wrapper.info["tree count"] == 2
    assert len(wrapper.trees) == 2

    mol_nodes = list(tree.molecules())
    rxn_nodes = list(tree.reactions())
    # Assert first tree
    assert wrapper.first_tree["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.first_tree["children"]) == 1

    child1 = wrapper.first_tree["children"][0]
    assert child1["smiles"] == rxn_nodes[0].smiles
    assert len(child1["children"]) == 2

    child_smiles = [child["smiles"] for child in child1["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles

    # Assert second tree
    assert wrapper.trees[1]["smiles"] == mol_nodes[0].smiles
    assert len(wrapper.trees[1]["children"]) == 1

    child1 = wrapper.trees[1]["children"][0]
    assert child1["smiles"] == rxn_nodes[0].smiles
    assert len(child1["children"]) == 2

    child_smiles = [child["smiles"] for child in child1["children"]]
    expected_smiles = [node.smiles for node in mol_nodes[1:]]
    assert child_smiles == expected_smiles[::-1]


def test_route_self_distance(load_reaction_tree):
    tree = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))
    wrapper = ReactionTreeWrapper(tree, exhaustive_limit=1)

    assert wrapper.distance_to(wrapper) == 0.0


def test_route_distances_random(load_reaction_tree):
    tree1 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))
    wrapper1 = ReactionTreeWrapper(tree1, exhaustive_limit=1)
    tree2 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 1))
    wrapper2 = ReactionTreeWrapper(tree2, exhaustive_limit=1)

    distances = list(wrapper1.distance_iter(wrapper2, exhaustive_limit=1))

    assert len(distances) == 2
    assert pytest.approx(distances[0], abs=1e-2) == 2.6522


def test_route_distances_exhaustive(load_reaction_tree):
    tree1 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))
    wrapper1 = ReactionTreeWrapper(tree1, exhaustive_limit=2)
    tree2 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 1))
    wrapper2 = ReactionTreeWrapper(tree2, exhaustive_limit=2)

    distances = list(wrapper1.distance_iter(wrapper2, exhaustive_limit=40))

    assert len(distances) == 2
    assert pytest.approx(distances[0], abs=1e-2) == 2.6522
    assert pytest.approx(min(distances), abs=1e-2) == 2.6522


def test_route_distances_semi_exhaustive(load_reaction_tree):
    tree1 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 0))
    wrapper1 = ReactionTreeWrapper(tree1, exhaustive_limit=1)
    tree2 = ReactionTree.from_dict(load_reaction_tree("routes_for_clustering.json", 1))
    wrapper2 = ReactionTreeWrapper(tree2, exhaustive_limit=2)

    distances = list(wrapper1.distance_iter(wrapper2, exhaustive_limit=1))

    assert len(distances) == 2
    assert pytest.approx(distances[0], abs=1e-2) == 2.6522
    assert pytest.approx(min(distances), abs=1e-2) == 2.6522
