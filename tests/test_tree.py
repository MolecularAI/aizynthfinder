import json

import pytest

from aizynthfinder.chem import TreeMolecule, Molecule, Reaction
from aizynthfinder.analysis import TreeAnalysis, ReactionTree, RouteCollection


@pytest.fixture
def setup_complete_tree(fresh_tree, mocker, mock_stock):
    tree = fresh_tree

    state1 = mocker.MagicMock()
    state1.mols = [
        TreeMolecule(
            parent=None,
            transform=0,
            smiles="CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1",
        )
    ]
    state1.in_stock_list = [False]
    state1.score = 0.049
    state1.is_solved = False
    node1 = mocker.MagicMock()
    node1.state = state1
    node1.parent = None
    node1.is_expanded = True
    action1 = (
        "([C:2]-[CH;D3;+0:1](-[C:3])-[C;H0;D3;+0:4](=[O;H0;D1;+0:6])-[c:5])"
        ">>(Cl-[CH;D3;+0:1](-[C:2])-[C:3]).(N#[C;H0;D2;+0:4]-[c:5]).([OH2;D0;+0:6])"
    )
    reaction1 = Reaction(state1.mols[0], action1)
    reaction1.apply()
    node1.__getitem__.return_value = {"action": reaction1}
    node1.is_terminal.return_value = False
    tree.root = node1

    state2 = mocker.MagicMock()
    state2.mols = [
        TreeMolecule(parent=state1.mols[0], smiles=smiles)
        for smiles in ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    ]
    state2.in_stock_list = [True, False, True]
    state2.score = 0.68
    state2.is_solved = False
    node2 = mocker.MagicMock()
    node2.parent = node1
    node2.is_expanded = True
    node2.state = state2
    node1.promising_child.return_value = node2
    node1.children.return_value = [node2]
    action2 = (
        "([O;D1;H0:2]=[C;H0;D3;+0:1](-[c:3])-[NH;D2;+0:4]-[c:5])"
        ">>(Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3]).([NH2;D1;+0:4]-[c:5])"
    )
    reaction2 = Reaction(state2.mols[1], action2)
    reaction2.apply()
    node2.__getitem__.return_value = {"action": reaction2}

    state3 = mocker.MagicMock()
    state3.mols = [
        TreeMolecule(parent=state2.mols[1], smiles=smiles)
        for smiles in ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    ]
    state3.in_stock_list = [True, True]
    state3.stock = mocker.MagicMock()
    state3.stock.__contains__.side_effect = [False, True, False, True, True, True]
    state3.score = 0.99
    state3.is_solved = True
    node3 = mocker.MagicMock()
    node3.parent = node2
    node3.state = state3
    node2.promising_child.return_value = node3
    node2.children.return_value = [node3]

    return tree, [node1, node2, node3]


def test_select_leaf_root(setup_complete_tree):
    tree, nodes = setup_complete_tree
    nodes[0].is_expanded = False

    leaf = tree.select_leaf()

    assert leaf is nodes[0]


def test_select_leaf(setup_complete_tree):
    tree, nodes = setup_complete_tree

    leaf = tree.select_leaf()

    assert leaf is nodes[2]


def test_backpropagation(setup_complete_tree):
    tree, nodes = setup_complete_tree
    score = 1.5

    tree.backpropagate(nodes[2], score)

    nodes[0].backpropagate.assert_called_once_with(nodes[1], score)
    nodes[1].backpropagate.assert_called_once_with(nodes[2], score)
    nodes[2].backpropagate.assert_not_called()


def test_route_to_node(setup_complete_tree):
    tree, nodes = setup_complete_tree

    routes, route_nodes = tree.route_to_node(nodes[2])

    assert len(routes) == 2
    assert len(nodes) == 3
    assert nodes[0] == route_nodes[0]
    assert nodes[1] == route_nodes[1]
    assert nodes[2] == route_nodes[2]


def test_create_graph(setup_complete_tree):
    tree, nodes = setup_complete_tree

    graph = tree.graph()

    assert len(graph) == 3
    assert list(graph.successors(nodes[0])) == [nodes[1]]
    assert list(graph.successors(nodes[1])) == [nodes[2]]


def test_sort_nodes(setup_complete_tree):
    tree, nodes = setup_complete_tree
    analysis = TreeAnalysis(tree)

    best_nodes = analysis.sort_nodes()

    assert len(best_nodes) == 3
    assert best_nodes[0].state.score == 0.99
    assert best_nodes[0] is nodes[2]

    best_nodes = analysis.sort_nodes(min_return=0)

    assert len(best_nodes) == 0


def test_route_to_reactiontree(setup_complete_tree):
    tree, nodes = setup_complete_tree
    analysis = TreeAnalysis(tree)

    reaction_tree = ReactionTree.from_analysis(analysis).graph

    reaction_nodes = [
        node.inchi_key for node in reaction_tree if isinstance(node, Molecule)
    ]
    tree_nodes = [
        mol.inchi_key
        for mol in nodes[0].state.mols + nodes[1].state.mols + nodes[2].state.mols
    ]
    assert len(reaction_nodes) == 6
    assert reaction_nodes == tree_nodes

    reaction_nodes = [node for node in reaction_tree if isinstance(node, Reaction)]
    assert len(reaction_nodes) == 2


def test_reactiontree_to_json(setup_complete_tree, shared_datadir):
    filename = str(shared_datadir / "sample_reaction.json")
    with open(filename, "r") as fileobj:
        expected = json.load(fileobj)
    tree, nodes = setup_complete_tree
    analysis = TreeAnalysis(tree)

    resp = ReactionTree.from_analysis(analysis).to_json()
    assert json.loads(resp) == expected


def test_reactiontree_from_dict(shared_datadir):
    filename = str(shared_datadir / "sample_reaction.json")
    with open(filename, "r") as fileobj:
        expected = json.load(fileobj)

    rt = ReactionTree.from_dict(expected)

    # Simply check that the to_dict() and from_dict() gives/produces the same dict
    resp = rt.to_dict()
    assert resp == expected


def test_create_route_collection(setup_complete_tree, mocker):
    tree, nodes = setup_complete_tree
    analysis = TreeAnalysis(tree)
    mocker.patch("aizynthfinder.analysis.ReactionTree.to_dict")
    mocker.patch("aizynthfinder.analysis.json.dumps")

    routes = RouteCollection.from_analysis(analysis, 5)

    assert len(routes) == 3
    assert routes[0]["score"] == 0.99
    assert routes[0]["node"] is nodes[2]
    reaction_nodes = [
        node for node in routes[0]["reaction_tree"].graph if isinstance(node, Reaction)
    ]
    assert len(reaction_nodes) == 2

    # Just see that the code does not crash, does not verify content
    assert len(routes.images) == 3
    assert len(routes.dicts) == 3
    assert len(routes.jsons) == 3
