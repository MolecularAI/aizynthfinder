from aizynthfinder.utils.serialization import MoleculeSerializer, MoleculeDeserializer
from aizynthfinder.chem import Molecule, TreeMolecule
from aizynthfinder.mcts.state import State
from aizynthfinder.mcts.node import Node
from aizynthfinder.mcts.mcts import SearchTree


def test_empty_store():
    serializer = MoleculeSerializer()

    assert serializer.store == {}


def test_add_single_mol():
    serializer = MoleculeSerializer()
    mol = Molecule(smiles="CCC")

    id_ = serializer[mol]

    assert id_ == id(mol)
    assert serializer.store == {id_: {"smiles": "CCC", "class": "Molecule"}}


def test_add_tree_mol():
    serializer = MoleculeSerializer()
    mol1 = TreeMolecule(parent=None, smiles="CCC", transform=1)
    mol2 = TreeMolecule(smiles="CCO", parent=mol1)

    id_ = serializer[mol2]

    assert id_ == id(mol2)
    assert list(serializer.store.keys()) == [id(mol1), id_]
    assert serializer.store == {
        id_: {
            "smiles": "CCO",
            "class": "TreeMolecule",
            "parent": id(mol1),
            "transform": 2,
        },
        id(mol1): {
            "smiles": "CCC",
            "class": "TreeMolecule",
            "parent": None,
            "transform": 1,
        },
    }


def test_deserialize_single_mol():
    store = {123: {"smiles": "CCC", "class": "Molecule"}}
    deserializer = MoleculeDeserializer(store)

    assert deserializer[123].smiles == "CCC"


def test_deserialize_tree_mols():
    store = {
        123: {
            "smiles": "CCC",
            "class": "TreeMolecule",
            "parent": None,
            "transform": 1,
        },
        234: {"smiles": "CCO", "class": "TreeMolecule", "parent": 123, "transform": 2},
    }

    deserializer = MoleculeDeserializer(store)

    assert deserializer[123].smiles == "CCC"
    assert deserializer[234].smiles == "CCO"
    assert deserializer[123].parent is None
    assert deserializer[234].parent is deserializer[123]
    assert deserializer[123].transform == 1
    assert deserializer[234].transform == 2


def test_chaining():
    serializer = MoleculeSerializer()
    mol1 = TreeMolecule(parent=None, smiles="CCC", transform=1)
    mol2 = TreeMolecule(smiles="CCO", parent=mol1)

    id_ = serializer[mol2]

    deserializer = MoleculeDeserializer(serializer.store)

    assert deserializer[id_].smiles == mol2.smiles
    assert deserializer[id(mol1)].smiles == mol1.smiles
    assert id(deserializer[id_]) != id_
    assert id(deserializer[id(mol1)]) != id(mol1)


def test_serialize_deserialize_state(default_config):
    mol = TreeMolecule(parent=None, smiles="CCC", transform=1)
    state0 = State([mol], default_config)
    serializer = MoleculeSerializer()

    state_serialized = state0.serialize(serializer)

    assert len(state_serialized["mols"]) == 1
    assert state_serialized["mols"][0] == id(mol)

    deserializer = MoleculeDeserializer(serializer.store)
    state1 = State.from_dict(state_serialized, default_config, deserializer)

    assert len(state1.mols) == 1
    assert state1.mols[0] == state0.mols[0]
    assert state1.in_stock_list == state0.in_stock_list
    assert state1.score == state0.score


def test_serialize_node(generate_root, simple_actions, mock_policy):
    serializer = MoleculeSerializer()
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, prior_list = mock_policy(root.state.mols[0])

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
    assert node_serialized["children_values"] == prior_list
    assert node_serialized["children_priors"] == prior_list
    assert node_serialized["children_visitations"] == [1, 1, 1]
    assert all(
        id(expected.mol) == actual["mol"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert all(
        expected.smarts == actual["smarts"]
        for expected, actual in zip(action_list, node_serialized["children_actions"])
    )
    assert all(
        expected.index == actual["index"]
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


def test_deserialize_node(generate_root, simple_actions, mock_policy, default_config):
    serializer = MoleculeSerializer()
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    action_list, prior_list = mock_policy(root.state.mols[0])
    root.expand()
    child = root.promising_child()
    node_serialized = root.serialize(serializer)
    deserializer = MoleculeDeserializer(serializer.store)

    root_new = Node.from_dict(node_serialized, None, default_config, deserializer)
    assert len(root_new.children()) == 1

    new_child = root_new.children()[0]
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
    fresh_tree,
    generate_root,
    simple_actions,
    mock_policy,
    default_config,
    mocker,
    tmpdir,
):
    serializer = MoleculeSerializer()
    root = generate_root("CCCCOc1ccc(CC(=O)N(C)O)cc1")
    fresh_tree.root = root
    action_list, prior_list = mock_policy(root.state.mols[0])
    root.expand()
    child = root.promising_child()
    mocked_json_dump = mocker.patch("aizynthfinder.mcts.mcts.json.dump")
    serializer = MoleculeSerializer()
    filename = str(tmpdir / "dummy.json")

    # Test serialization

    fresh_tree.serialize(filename)

    expected_dict = {"tree": root.serialize(serializer), "molecules": serializer.store}
    mocked_json_dump.assert_called_once_with(
        expected_dict, mocker.ANY, indent=mocker.ANY
    )

    # Test deserialization

    mocker.patch("aizynthfinder.mcts.mcts.json.load", return_value=expected_dict)

    new_tree = SearchTree.from_json(filename, default_config)
    root_new = new_tree.root
    assert len(root_new.children()) == 1

    new_child = root_new.children()[0]
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
