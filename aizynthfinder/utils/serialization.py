""" Module containing helper classes and routines for serialization.
"""
import aizynthfinder.chem


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

    def __init__(self):
        self._store = {}

    def __getitem__(self, mol):
        if mol is None:
            return None

        id_ = id(mol)
        if id_ not in self._store:
            self._add_mol(mol)
        return id_

    @property
    def store(self):
        """ Return all serialized molecules as a dictionary
        """
        return self._store

    def _add_mol(self, mol):
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

    def __init__(self, store):
        self._objects = {}
        self._create_molecules(store)

    def __getitem__(self, id_):
        if id_ is None:
            return None
        return self._objects[id_]

    def _create_molecules(self, store):
        for id_, spec in store.items():
            if isinstance(id_, str):
                id_ = int(id_)

            cls = spec["class"]
            if "parent" in spec:
                spec["parent"] = self[spec["parent"]]

            kwargs = dict(spec)
            del kwargs["class"]
            self._objects[id_] = getattr(aizynthfinder.chem, cls)(**kwargs)
