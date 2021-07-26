from aizynthfinder.chem.serialization import MoleculeSerializer, MoleculeDeserializer
from aizynthfinder.chem import TreeMolecule
from aizynthfinder.search.mcts import MctsState
from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.search.mcts import MctsSearchTree


def test_serialize_deserialize_state(default_config):
    mol = TreeMolecule(parent=None, smiles="CCC", transform=1)
    state0 = MctsState([mol], default_config)
    serializer = MoleculeSerializer()

    state_serialized = state0.serialize(serializer)

    assert len(state_serialized["mols"]) == 1
    assert state_serialized["mols"][0] == id(mol)

    deserializer = MoleculeDeserializer(serializer.store)
    state1 = MctsState.from_dict(state_serialized, default_config, deserializer)

    assert len(state1.mols) == 1
    assert state1.mols[0] == state0.mols[0]
    assert state1.in_stock_list == state0.in_stock_list
    assert state1.score == state0.score


def test_serialize_node(setup_mcts_search):
    serializer = MoleculeSerializer()
    root, _, _ = setup_mcts_search

    state_serialized = root.state.serialize(serializer)
    node_serialized = root.serialize(serializer)
    assert not node_serialized["is_expanded"]
    assert node_serialized["state"] == state_serialized
    assert node_serialized["children_values"] == []
    assert node_serialized["children_priors"] == []
    assert node_serialized["children_visitations"] == []
    assert node_serialized["children"] == []

    root.expand()

    node_serialized = root.serialize(serializer)
    action_list = root.children_view()["actions"]
    assert node_serialized["children_values"] == [0.7, 0.5, 0.3]
    assert node_serialized["children_priors"] == [0.7, 0.5, 0.3]
    assert node_serialized["children_visitations"] == [1, 1, 1]
    assert all(
        id(expected.mol) == actual["mol"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert all(
        expected.reactants_str == actual["reactants_str"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert all(
        expected.index == actual["index"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert all(
        expected.metadata == actual["metadata"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert node_serialized["children"] == [None, None, None]
    assert node_serialized["is_expanded"]

    child = root.promising_child()

    node_serialized = root.serialize(serializer)
    state_serialized = child.state.serialize(serializer)
    assert node_serialized["is_expanded"]
    assert node_serialized["children"][1] is None
    assert node_serialized["children"][2] is None
    assert node_serialized["children"][0]["state"] == state_serialized
    assert node_serialized["children"][0]["children_values"] == []
    assert node_serialized["children"][0]["children_priors"] == []
    assert node_serialized["children"][0]["children_visitations"] == []
    assert node_serialized["children"][0]["children"] == []
    assert not node_serialized["children"][0]["is_expanded"]


def test_deserialize_node(setup_mcts_search, default_config):
    serializer = MoleculeSerializer()
    root, _, _ = setup_mcts_search
    root.expand()
    child = root.promising_child()
    node_serialized = root.serialize(serializer)
    deserializer = MoleculeDeserializer(serializer.store)

    root_new = MctsNode.from_dict(node_serialized, None, default_config, deserializer)
    assert len(root_new.children) == 1

    new_child = root_new.children[0]
    assert root_new.children_view()["values"] == root.children_view()["values"]
    assert root_new.children_view()["priors"] == root.children_view()["priors"]
    assert (
        root_new.children_view()["visitations"] == root.children_view()["visitations"]
    )
    assert root_new.is_expanded
    assert new_child.children_view()["values"] == child.children_view()["values"]
    assert new_child.children_view()["priors"] == child.children_view()["priors"]
    assert (
        new_child.children_view()["visitations"] == child.children_view()["visitations"]
    )
    assert not new_child.is_expanded
    assert str(root_new.state) == str(root.state)
    assert str(new_child.state) == str(child.state)


def test_serialize_deserialize_tree(
    setup_complete_mcts_tree,
    default_config,
    mocker,
    tmpdir,
):
    tree, nodes = setup_complete_mcts_tree
    root, child, _ = nodes
    mocked_json_dump = mocker.patch("aizynthfinder.search.mcts.search.json.dump")
    serializer = MoleculeSerializer()
    filename = str(tmpdir / "dummy.json")

    # Test serialization

    tree.serialize(filename)

    expected_dict = {"tree": root.serialize(serializer), "molecules": serializer.store}
    mocked_json_dump.assert_called_once_with(
        expected_dict, mocker.ANY, indent=mocker.ANY
    )

    # Test deserialization

    mocker.patch(
        "aizynthfinder.search.mcts.search.json.load", return_value=expected_dict
    )

    new_tree = MctsSearchTree.from_json(filename, default_config)
    root_new = new_tree.root
    assert len(root_new.children) == 1

    new_child = root_new.children[0]
    assert root_new.children_view()["values"] == root.children_view()["values"]
    assert root_new.children_view()["priors"] == root.children_view()["priors"]
    assert (
        root_new.children_view()["visitations"] == root.children_view()["visitations"]
    )
    assert root_new.is_expanded
    assert new_child.children_view()["values"] == child.children_view()["values"]
    assert new_child.children_view()["priors"] == child.children_view()["priors"]
    assert (
        new_child.children_view()["visitations"] == child.children_view()["visitations"]
    )
    assert new_child.is_expanded
    assert str(root_new.state) == str(root.state)
    assert str(new_child.state) == str(child.state)
