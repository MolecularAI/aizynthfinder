import pytest

from aizynthfinder.reactiontree import ReactionTree


def test_mcts_route_to_reactiontree(setup_linear_mcts, load_reaction_tree):
    def remove_metadata(tree_dict):
        if "metadata" in tree_dict:
            tree_dict["metadata"] = {}
        for child in tree_dict.get("children", []):
            remove_metadata(child)

    _, node = setup_linear_mcts()
    expected_dict = load_reaction_tree("linear_route.json")

    reaction_tree = node.to_reaction_tree()
    assert reaction_tree.metadata == {"created_at_iteration": 1, "is_solved": True}

    mol_nodes = list(reaction_tree.molecules())
    assert len(mol_nodes) == 5

    reaction_nodes = list(reaction_tree.reactions())
    assert len(reaction_nodes) == 2

    leaf_nodes = list(reaction_tree.leafs())
    assert len(leaf_nodes) == 3
    assert set(node.smiles for node in leaf_nodes) == {
        "NC1CCCC(C2C=CC=C2)C1",
        "c1ccccc1",
        "OOc1ccccc1",
    }
    actual_dict = reaction_tree.to_dict(include_metadata=True)
    remove_metadata(actual_dict)
    assert actual_dict == expected_dict


def test_reactiontree_from_dict(load_reaction_tree):
    expected = load_reaction_tree("linear_route.json")

    rt = ReactionTree.from_dict(expected)

    # Simply check that the to_dict() and from_dict() gives/produces the same dict
    resp = rt.to_dict(include_metadata=True)
    assert resp == expected


def test_reactiontree_to_image(load_reaction_tree, mocker):
    patched_make_image = mocker.patch("aizynthfinder.reactiontree.RouteImageFactory")

    tree = load_reaction_tree("linear_route.json")
    rt = ReactionTree.from_dict(tree)

    rt.to_image()

    patched_make_image.assert_called_once()


def test_route_node_depth_from_mcts(setup_branched_reaction_tree):
    rt = setup_branched_reaction_tree()

    mols = list(rt.molecules())
    assert rt.depth(mols[0]) == 0
    assert rt.depth(mols[1]) == 2
    assert rt.depth(mols[2]) == 2
    assert rt.depth(mols[3]) == 4
    assert rt.depth(mols[4]) == 4
    assert rt.depth(mols[5]) == 6
    assert rt.depth(mols[6]) == 6
    assert rt.depth(mols[7]) == 4
    assert rt.depth(mols[8]) == 4

    rxns = list(rt.reactions())
    assert rt.depth(rxns[0]) == 1
    assert rt.depth(rxns[1]) == 3
    assert rt.depth(rxns[2]) == 5
    assert rt.depth(rxns[3]) == 3

    for mol in rt.molecules():
        assert rt.depth(mol) == 2 * rt.graph.nodes[mol]["transform"]


def test_route_node_depth_from_json(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json")
    rt = ReactionTree.from_dict(dict_)

    # Molecules loaded in a different order than when created from MCTS

    mols = list(rt.molecules())
    assert rt.depth(mols[0]) == 0
    assert rt.depth(mols[1]) == 2
    assert rt.depth(mols[2]) == 4
    assert rt.depth(mols[3]) == 4
    assert rt.depth(mols[4]) == 6
    assert rt.depth(mols[5]) == 6
    assert rt.depth(mols[6]) == 2
    assert rt.depth(mols[7]) == 4
    assert rt.depth(mols[8]) == 4

    rxns = list(rt.reactions())
    assert rt.depth(rxns[0]) == 1
    assert rt.depth(rxns[1]) == 3
    assert rt.depth(rxns[2]) == 5
    assert rt.depth(rxns[3]) == 3

    for mol in rt.molecules():
        assert rt.depth(mol) == 2 * rt.graph.nodes[mol]["transform"]


def test_route_distance_self(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json", 0, remove_metadata=False)
    rt = ReactionTree.from_dict(dict_)

    assert rt.distance_to(rt) == 0.0


def test_route_distance_other(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json", remove_metadata=False)
    rt1 = ReactionTree.from_dict(dict_)
    dict_ = load_reaction_tree("linear_route.json", remove_metadata=False)
    rt2 = ReactionTree.from_dict(dict_)

    dist = rt1.distance_to(rt2)

    assert pytest.approx(dist, abs=1e-2) == 0.6028


@pytest.mark.parametrize(
    ("filename", "expected"),
    [("linear_route.json", False), ("branched_route.json", True)],
)
def test_route_is_branched(load_reaction_tree, filename, expected):
    rt = ReactionTree.from_dict(load_reaction_tree(filename))

    assert rt.is_branched() == expected


def test_route_hash(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json")
    rt = ReactionTree.from_dict(dict_)

    assert rt.hash_key() == "1514c305b6dcddc0a2e4133ff77cd53893d25c4e5caaca7dc53490fe"


def test_subtrees(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json")
    rt = ReactionTree.from_dict(dict_)

    subtrees = list(rt.subtrees())

    assert len(subtrees) == 3

    mols = list(subtrees[0].molecules())
    assert len(mols) == 5
    assert [subtrees[0].depth(mol) for mol in mols] == [0, 2, 2, 4, 4]

    mols = list(subtrees[1].molecules())
    assert len(mols) == 3
    assert [subtrees[1].depth(mol) for mol in mols] == [0, 2, 2]

    mols = list(subtrees[2].molecules())
    assert len(mols) == 3
    assert [subtrees[2].depth(mol) for mol in mols] == [0, 2, 2]


def test_reactiontree_parent_mol(load_reaction_tree):
    dict_ = load_reaction_tree("linear_route.json")
    rt = ReactionTree.from_dict(dict_)
    molecules = list(rt.molecules())

    assert rt.parent_molecule(molecules[1]) is molecules[0]
    assert rt.parent_molecule(molecules[2]) is molecules[0]

    with pytest.raises(ValueError):
        rt.parent_molecule(molecules[0])


def test_reactiontree_child_reactions(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json")
    rt = ReactionTree.from_dict(dict_)
    reactions = list(rt.reactions())

    child_reactions = rt.child_reactions(reactions[0])
    child_reactions2 = rt.child_reactions(reactions[1])

    assert len(child_reactions) == 2
    assert len(child_reactions2) == 1

    assert child_reactions[0] is reactions[1]
    assert child_reactions[1] is reactions[3]

    assert child_reactions2[0] is reactions[2]
    assert len(rt.child_reactions(reactions[-1])) == 0
