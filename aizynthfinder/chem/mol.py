""" Module containing classes to deal with Molecules - mostly wrappers around rdkit routines.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

from aizynthfinder.utils.exceptions import MoleculeException


if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Dict,
        Optional,
        Union,
        Tuple,
        RdMol,
        Sequence,
        List,
        Callable,
    )


class Molecule:
    """
    A base class for molecules. Encapsulate an RDKit mol object and
    functions that can be applied to such a molecule.

    The objects of this class is hashable by the inchi key and hence
    comparable with the equality operator.

    :ivar rd_mol: the RDkit mol object that is encapsulated
    :ivar smiles: the SMILES representation of the molecule

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(
        self, rd_mol: RdMol = None, smiles: str = None, sanitize: bool = False
    ) -> None:
        if not rd_mol and not smiles:
            raise MoleculeException(
                "Need to provide either a rdkit Mol object or smiles to create a molecule"
            )

        if rd_mol:
            self.rd_mol = rd_mol
            self.smiles = Chem.MolToSmiles(rd_mol)
        else:
            self.smiles = smiles
            self.rd_mol = Chem.MolFromSmiles(smiles, sanitize=False)

        self._inchi_key: Optional[str] = None
        self._inchi: Optional[str] = None
        self._fingerprints: Dict[Union[Tuple[int, int], Tuple[int]], np.ndarray] = {}
        self._is_sanitized: bool = False

        # Atom mapping -> atom index dictionary
        self._atom_mappings: Dict[int, int] = {}
        # Atom index -> atom mapping dictionary
        self._reverse_atom_mappings: Dict[int, int] = {}

        if sanitize:
            self.sanitize()

    def __hash__(self) -> int:
        return hash(self.inchi_key)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Molecule):
            return False
        return self.inchi_key == other.inchi_key

    def __len__(self) -> int:
        return self.rd_mol.GetNumAtoms()

    def __str__(self) -> str:
        return self.smiles

    @property
    def inchi(self) -> str:
        """
        The inchi representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi
        """
        if not self._inchi:
            self.sanitize(raise_exception=False)
            self._inchi = Chem.MolToInchi(self.rd_mol)
            if self._inchi is None:
                raise MoleculeException("Could not make InChI")
        return self._inchi

    @property
    def inchi_key(self) -> str:
        """
        The inchi key representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi key
        """
        if not self._inchi_key:
            self.sanitize(raise_exception=False)
            self._inchi_key = Chem.MolToInchiKey(self.rd_mol)
            if self._inchi_key is None:
                raise MoleculeException("Could not make InChI key")
        return self._inchi_key

    @property
    def index_to_mapping(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom indices to atom mappings"""
        if not self._reverse_atom_mappings:
            self._reverse_atom_mappings = {
                index: mapping for mapping, index in self.mapping_to_index.items()
            }
        return self._reverse_atom_mappings

    @property
    def mapping_to_index(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom mappings to atom indices"""
        if not self._atom_mappings:
            self._atom_mappings = {
                atom.GetAtomMapNum(): atom.GetIdx()
                for atom in self.rd_mol.GetAtoms()
                if atom.GetAtomMapNum()
            }
        return self._atom_mappings

    @property
    def weight(self) -> float:
        """Return the exact molecular weight of the molecule"""
        self.sanitize(raise_exception=False)
        return Descriptors.ExactMolWt(self.rd_mol)

    def basic_compare(self, other: "Molecule") -> bool:
        """
        Compare this molecule to another but only to
        the basic part of the inchi key, thereby ignoring stereochemistry etc

        :param other: the molecule to compare to
        :return: True if chemical formula and connectivity is the same
        """
        return self.inchi_key[:14] == other.inchi_key[:14]

    def fingerprint(self, radius: int, nbits: int = 2048) -> np.ndarray:
        """
        Returns the Morgan fingerprint of the molecule

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint
        :return: the fingerprint
        """
        key = radius, nbits

        if key not in self._fingerprints:
            self.sanitize()
            bitvect = AllChem.GetMorganFingerprintAsBitVect(self.rd_mol, *key)
            array = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(bitvect, array)
            self._fingerprints[key] = array

        return self._fingerprints[key]

    def has_atom_mapping(self) -> bool:
        """
        Determines if a the molecule has atom mappings

        :return: True if at least one atom has a mapping
        """
        for atom in self.rd_mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                return True
        return False

    def make_unique(self) -> "UniqueMolecule":
        """
        Returns an instance of the UniqueMolecule class that
        is representing the same molecule but is not hashable or comparable.

        :return: the unique molecule
        """
        return UniqueMolecule(rd_mol=self.rd_mol)

    def remove_atom_mapping(self, exceptions: Sequence[int] = None) -> None:
        """
        Remove all mappings of the atoms and update the smiles

        :param exceptions: keep the listed atom mappings
        """
        exceptions = exceptions or []
        for atom in self.rd_mol.GetAtoms():
            if exceptions and atom.GetAtomMapNum() in exceptions:
                continue
            atom.SetAtomMapNum(0)
        self.smiles = Chem.MolToSmiles(self.rd_mol)
        self._clear_cache()

    def sanitize(self, raise_exception: bool = True) -> None:
        """
        Sanitizes the molecule if it has not been done before.

        :param raise_exception: if True will raise exception on failed sanitation
        :raises MoleculeException: if the molecule could not be sanitized
        """
        if self._is_sanitized:
            return

        try:
            AllChem.SanitizeMol(self.rd_mol)
        # pylint: disable=bare-except
        except:  # noqa, there could be many reasons why the molecule cannot be sanitized
            if raise_exception:
                raise MoleculeException(f"Unable to sanitize molecule ({self.smiles})")
            self.rd_mol = Chem.MolFromSmiles(self.smiles, sanitize=False)
            return

        self.smiles = Chem.MolToSmiles(self.rd_mol)
        self._clear_cache()
        self._is_sanitized = True

    def _clear_cache(self):
        self._inchi = None
        self._inchi_key = None
        self._fingerprints = {}
        self._atom_mappings = {}
        self._reverse_atom_mappings = {}


class TreeMolecule(Molecule):
    """
    A special molecule that keeps a reference to a parent molecule.

    If the class is instantiated without specifying the `transform` argument,
    it is computed by increasing the value of the `parent.transform` variable.

    If no parent is provided the atoms with atom mapping number are tracked
    and inherited to children.

    :ivar mapped_mol: the tracked molecule with atom mappings
    :ivar mapped_smiles: the SMILES of the tracked molecule with atom mappings
    :ivar original_smiles: the SMILES as passed when instantiating the class
    :ivar parent: parent molecule
    :ivar transform: a numerical number corresponding to the depth in the tree
    :ivar tracked_atom_indices: tracked atom indices and what indices they correspond to in this molecule

    :param parent: a TreeMolecule object that is the parent
    :param transform: the transform value, defaults to None
    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :param mapping_update_callback: if given will call this method before setting up the `mapped_smiles`
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    # pylint: disable=too-many-arguments
    def __init__(
        self,
        parent: Optional["TreeMolecule"],
        transform: int = None,
        rd_mol: RdMol = None,
        smiles: str = None,
        sanitize: bool = False,
        mapping_update_callback: Callable[["TreeMolecule"], None] = None,
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)
        self.parent = parent
        if transform is None and parent and parent.transform is not None:
            self.transform: int = parent.transform + 1
        else:
            self.transform = transform or 0

        self.original_smiles = smiles
        self.tracked_atom_indices: Dict[int, Optional[int]] = {}
        self.mapped_mol = Chem.Mol(self.rd_mol)
        if not self.parent:
            self._init_tracking()
        elif mapping_update_callback is not None:
            mapping_update_callback(self)

        AllChem.SanitizeMol(self.mapped_mol)
        self.mapped_smiles = Chem.MolToSmiles(self.mapped_mol)

        if self.parent:
            self.remove_atom_mapping()
            self._update_tracked_atoms()

    @property
    def mapping_to_index(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom mappings to atom indices"""
        if not self._atom_mappings:
            self._atom_mappings = {
                atom.GetAtomMapNum(): atom.GetIdx()
                for atom in self.mapped_mol.GetAtoms()
                if atom.GetAtomMapNum()
            }
        return self._atom_mappings

    def _init_tracking(self):
        self.tracked_atom_indices = dict(self.mapping_to_index)
        for idx, atom in enumerate(self.mapped_mol.GetAtoms()):
            atom.SetAtomMapNum(idx + 1)
        self._atom_mappings = {}

    def _update_tracked_atoms(self) -> None:
        if self.parent is None:
            return

        if not self.parent.tracked_atom_indices:
            return

        parent2child_map = {
            atom_index: self.mapping_to_index.get(mapping_index)
            for mapping_index, atom_index in self.parent.mapping_to_index.items()
        }

        self.tracked_atom_indices = {
            tracked_index: parent2child_map[parent_index]  # type: ignore
            for tracked_index, parent_index in self.parent.tracked_atom_indices.items()
        }


class UniqueMolecule(Molecule):
    """
    A special molecule with the hash set to the id of the object.
    Therefore no two instances of this class will be comparable.

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(
        self, rd_mol: RdMol = None, smiles: str = None, sanitize: bool = False
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, _) -> bool:
        return False


def none_molecule() -> UniqueMolecule:
    """Return an empty molecule"""
    return UniqueMolecule(rd_mol=Chem.MolFromSmiles(""))
