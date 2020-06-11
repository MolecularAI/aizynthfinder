def test_root_state_properties(generate_root):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    root2 = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")

    assert round(root.state.score, 4) == 0.0491
    assert root.state.stock_availability == ["Not in stock"]
    assert hash(root.state) == hash(root2.state)


def test_expand_root_node(generate_root, mock_policy):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, prior_list = mock_policy(root.state.mols[0])

    root.expand()

    view = root.children_view()
    assert view["actions"] == action_list
    assert view["priors"] == prior_list
    assert view["values"] == prior_list
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_expand_root_with_default_priors(generate_root, set_default_prior, mock_policy):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, prior_list = mock_policy(root.state.mols[0])
    set_default_prior(0.01)

    root.expand()

    view = root.children_view()
    assert view["actions"] == action_list
    assert view["priors"] == prior_list
    assert view["values"] == [0.01, 0.01, 0.01]
    assert view["visitations"] == [1, 1, 1]
    assert view["objects"] == [None, None, None]


def test_expand_when_solved(generate_root, mock_policy, mock_stock):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    mock_policy(root.state.mols[0])
    root.expand()
    child = root.promising_child()
    mock_stock([True, True])

    child.expand()

    assert not hasattr(child, "children_actions")
    assert child.is_terminal


def test_promising_child_of_root(generate_root, simple_actions, mock_policy):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, _ = mock_policy(root.state.mols[0])
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][0] is child
    assert root.children() == [child]
    assert root[child]["visitations"] == 1
    assert root[child]["value"] == 0.7
    assert root[child]["action"] is action_list[0]
    assert root[child]["prior"] == 0.7


def test_promising_child_of_root_invalid_action(generate_root, simple_actions, mocker):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    mocked_get_action = mocker.patch("aizynthfinder.mcts.policy.Policy.get_actions")
    action_list, prior_list = simple_actions(root.state.mols[0])
    prior_list[2] = 0.9  # This will select the invalid action
    mocked_get_action.return_value = action_list, prior_list
    root.expand()

    child = root.promising_child()

    view = root.children_view()
    assert view["objects"][0] is child
    assert child is not None
    assert root.children() == [child]


def test_backpropagate(generate_root, simple_actions, mock_policy):
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, _ = mock_policy(root.state.mols[0])
    root.expand()
    child = root.promising_child()
    view_prior = root.children_view()

    root.backpropagate(child, 1.5)

    view_post = root.children_view()
    assert view_post["visitations"][0] == view_prior["visitations"][0] + 1
    assert view_prior["visitations"][1:] == view_post["visitations"][1:]
    assert view_post["values"][0] == view_prior["values"][0] + 1.5
    assert view_prior["values"][1:] == view_post["values"][1:]
