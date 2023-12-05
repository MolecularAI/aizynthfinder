from aizynthfinder.chem.serialization import MoleculeSerializer, MoleculeDeserializer
from aizynthfinder.chem import Molecule, TreeMolecule


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
