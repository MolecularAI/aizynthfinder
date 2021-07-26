""" Module containing helper classes and routines for serialization.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import aizynthfinder.chem
from aizynthfinder.utils.loading import load_dynamic_class

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Optional, Sequence, Dict, Any, StrDict
    from aizynthfinder.chem import RetroReaction


class MoleculeSerializer:
    """
    Utility class for serializing molecules

    The id of the molecule to be serialized can be obtained with:

    .. code-block::

        serializer = MoleculeSerializer()
        mol = Molecule(smiles="CCCO")
        idx = serializer[mol]

    which will take care of the serialization of the molecule.
    """

    def __init__(self) -> None:
        self._store: Dict[int, Any] = {}

    def __getitem__(self, mol: Optional[aizynthfinder.chem.Molecule]) -> Optional[int]:
        if mol is None:
            return None

        id_ = id(mol)
        if id_ not in self._store:
            self._add_mol(mol)
        return id_

    @property
    def store(self) -> Dict[int, Any]:
        """Return all serialized molecules as a dictionary"""
        return self._store

    def _add_mol(self, mol: aizynthfinder.chem.Molecule) -> None:
        id_ = id(mol)
        dict_ = {"smiles": mol.smiles, "class": mol.__class__.__name__}
        if isinstance(mol, aizynthfinder.chem.TreeMolecule):
            dict_["parent"] = self[mol.parent]
            dict_["transform"] = mol.transform
        self._store[id_] = dict_


class MoleculeDeserializer:
    """
    Utility class for deserializing molecules.
    The serialized molecules are created upon instantiation of the class.

    The deserialized molecules can be obtained with:

    .. code-block::

        deserializer = MoleculeDeserializer()
        mol = deserializer[idx]

    """

    def __init__(self, store: Dict[int, Any]) -> None:
        self._objects: Dict[int, Any] = {}
        self._create_molecules(store)

    def __getitem__(self, id_: Optional[int]) -> Optional[aizynthfinder.chem.Molecule]:
        if id_ is None:
            return None
        return self._objects[id_]

    def get_tree_molecules(
        self, ids: Sequence[int]
    ) -> Sequence[aizynthfinder.chem.TreeMolecule]:
        """
        Return multiple deserialized tree molecules

        :param ids: the list of IDs to deserialize
        :return: the molecule objects
        """
        objects = []
        for id_ in ids:
            obj = self[id_]
            if obj is None or not isinstance(obj, aizynthfinder.chem.TreeMolecule):
                raise ValueError(f"Failed to deserialize molecule with id {id_}")
            objects.append(obj)
        return objects

    def _create_molecules(self, store: dict) -> None:
        for id_, spec in store.items():
            if isinstance(id_, str):
                id_ = int(id_)

            cls = spec["class"]
            if "parent" in spec:
                spec["parent"] = self[spec["parent"]]

            kwargs = dict(spec)
            del kwargs["class"]
            self._objects[id_] = getattr(aizynthfinder.chem, cls)(**kwargs)


def serialize_action(
    action: RetroReaction, molecule_store: MoleculeSerializer
) -> StrDict:
    """
    Serialize a retrosynthesis action

    :param action: the (re)action to serialize
    :param molecule_store: the molecule serialization object
    :return: the action as a dictionary
    """
    dict_ = action.to_dict()
    dict_["mol"] = molecule_store[dict_["mol"]]
    dict_["class"] = f"{action.__class__.__module__}.{action.__class__.__name__}"
    return dict_


def deserialize_action(
    dict_: StrDict, molecule_store: MoleculeDeserializer
) -> RetroReaction:
    """
    Deserialize a retrosynthesis action

    :param dict_: the (re)action as a dictionary
    :param molecule_store: the molecule deserialization object
    :return: the created action object
    """
    mol_spec = dict_.pop("mol")
    mol = molecule_store.get_tree_molecules([mol_spec])[0]
    try:
        class_spec = dict_.pop("class")
    except KeyError:
        class_spec = "aizynthfinder.chem.TemplatedRetroReaction"
    cls = load_dynamic_class(class_spec)
    return cls(mol, **dict_)
