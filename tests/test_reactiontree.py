import pytest

from aizynthfinder.reactiontree import ReactionTree


def test_mcts_route_to_reactiontree(setup_linear_mcts, load_reaction_tree):
    _, node = setup_linear_mcts()
    expected_dict = load_reaction_tree("linear_route.json")

    reaction_tree = node.to_reaction_tree()

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
    assert reaction_tree.to_dict() == expected_dict


def test_reactiontree_from_dict(load_reaction_tree):
    expected = load_reaction_tree("linear_route.json")

    rt = ReactionTree.from_dict(expected)

    # Simply check that the to_dict() and from_dict() gives/produces the same dict
    resp = rt.to_dict()
    assert resp == expected


def test_reactiontree_to_image(load_reaction_tree, mocker):
    patched_make_image = mocker.patch("aizynthfinder.reactiontree.make_graphviz_image")

    tree = load_reaction_tree("linear_route.json")
    rt = ReactionTree.from_dict(tree)

    rt.to_image()

    patched_make_image.assert_called_once()
    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules()))
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions()))


def test_reactiontree_to_image_hiding(load_reaction_tree, mocker):
    patched_make_image = mocker.patch("aizynthfinder.reactiontree.make_graphviz_image")

    tree = load_reaction_tree("linear_route.json", 1)
    rt = ReactionTree.from_dict(tree)
    for idx, reaction in enumerate(rt.reactions()):
        if idx == 0:
            continue
        rt.graph.nodes[reaction]["hide"] = True
    for idx, mol in enumerate(rt.leafs()):
        if idx == 2:
            continue
        rt.graph.nodes[mol]["hide"] = True

    rt.to_image(show_all=True)

    patched_make_image.assert_called_once()
    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules()))
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions()))

    patched_make_image.reset_mock()

    rt.to_image(show_all=False)

    assert len(patched_make_image.call_args[0][0]) == len(list(rt.molecules())) - 2
    assert len(patched_make_image.call_args[0][1]) == len(list(rt.reactions())) - 1


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
    dict_ = load_reaction_tree("branched_route.json", 0)
    rt = ReactionTree.from_dict(dict_)

    assert rt.distance_to(rt) == 0.0


def test_route_distance_other(load_reaction_tree):
    dict_ = load_reaction_tree("branched_route.json")
    rt1 = ReactionTree.from_dict(dict_)
    dict_ = load_reaction_tree("linear_route.json")
    rt2 = ReactionTree.from_dict(dict_)

    dist = rt1.distance_to(rt2, content="molecules")

    assert pytest.approx(dist, abs=1e-2) == 4.000


@pytest.mark.parametrize(
    ("filename", "expected"),
    [("linear_route.json", False), ("branched_route.json", True)],
)
def test_route_is_branched(load_reaction_tree, filename, expected):
    rt = ReactionTree.from_dict(load_reaction_tree(filename))

    assert rt.is_branched() == expected
