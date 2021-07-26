def test_select_leaf_root(setup_complete_mcts_tree):
    tree, nodes = setup_complete_mcts_tree
    nodes[0].is_expanded = False

    leaf = tree.select_leaf()

    assert leaf is nodes[0]


def test_select_leaf(setup_complete_mcts_tree):
    tree, nodes = setup_complete_mcts_tree

    leaf = tree.select_leaf()

    assert leaf is nodes[2]


def test_backpropagation(setup_complete_mcts_tree, mocker):
    tree, nodes = setup_complete_mcts_tree
    for node in nodes:
        node.backpropagate = mocker.MagicMock()
    score = 1.5

    tree.backpropagate(nodes[2], score)

    nodes[0].backpropagate.assert_called_once_with(nodes[1], score)
    nodes[1].backpropagate.assert_called_once_with(nodes[2], score)
    nodes[2].backpropagate.assert_not_called()


def test_route_to_node(setup_complete_mcts_tree):
    tree, nodes = setup_complete_mcts_tree

    actions, route_nodes = nodes[2].path_to()

    assert len(actions) == 2
    assert len(nodes) == 3
    assert nodes[0] == route_nodes[0]
    assert nodes[1] == route_nodes[1]
    assert nodes[2] == route_nodes[2]


def test_create_graph(setup_complete_mcts_tree):
    tree, nodes = setup_complete_mcts_tree

    graph = tree.graph()

    assert len(graph) == 3
    assert list(graph.successors(nodes[0])) == [nodes[1]]
    assert list(graph.successors(nodes[1])) == [nodes[2]]
