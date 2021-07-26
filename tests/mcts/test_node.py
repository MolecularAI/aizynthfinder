def test_root_state_properties(generate_root):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    root2 = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")

    assert round(root.state.score, 4) == 0.0491
    assert root.state.stock_availability == ["Not in stock"]
    assert hash(root.state) == hash(root2.state)


def test_expand_root_node(setup_mcts_search):
    root, _, _ = setup_mcts_search

    root.expand()

    view = root.children_view()
    assert len(view["actions"]) == 3
    assert view["priors"] == [0.7, 0.5, 0.3]
    assert view["values"] == [0.7, 0.5, 0.3]
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_expand_root_with_default_priors(setup_mcts_search, set_default_prior):
    root, _, _ = setup_mcts_search
    set_default_prior(0.01)

    root.expand()

    view = root.children_view()
    assert len(view["actions"]) == 3
    assert view["priors"] == [0.7, 0.5, 0.3]
    assert view["values"] == [0.01, 0.01, 0.01]
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_expand_when_solved(setup_mcts_search, setup_stock):
    root, _, _ = setup_mcts_search
    root.expand()
    child = root.promising_child()
    setup_stock(None, root.state.mols[0])

    child.expand()

    assert not hasattr(child, "children_actions")
    assert child.is_terminal


def test_promising_child_of_root(setup_mcts_search):
    root, _, _ = setup_mcts_search
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][0] is child
    assert root.children == [child]
    assert root[child]["visitations"] == 1
    assert root[child]["value"] == 0.7
    assert root[child]["prior"] == 0.7


def test_promising_child_of_root_invalid_action(setup_mcts_search):
    root, strategy, _ = setup_mcts_search
    strategy.lookup[root.state.mols[0].smiles][2]["prior"] = 0.9
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][0] is child
    assert child is not None
    assert root.children == [child]


def test_promising_child_with_filter(setup_mcts_search):
    root, _, strategy = setup_mcts_search
    strategy.lookup["CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"] = 0.8
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][0] is child
    assert child is not None
    assert root.children == [child]


def test_promising_child_with_filter_reject(setup_mcts_search):
    root, _, strategy = setup_mcts_search
    strategy.lookup["CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"] = 0.0
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][1] is child
    assert view["values"] == [-1000000.0, 0.5, 0.3]


def test_backpropagate(setup_mcts_search):
    root, _, _ = setup_mcts_search
    root.expand()
    child = root.promising_child()
    view_prior = root.children_view()

    root.backpropagate(child, 1.5)

    view_post = root.children_view()
    assert view_post["visitations"][0] == view_prior["visitations"][0] + 1
    assert view_prior["visitations"][1:] == view_post["visitations"][1:]
    assert view_post["values"][0] == view_prior["values"][0] + 1.5
    assert view_prior["values"][1:] == view_post["values"][1:]
