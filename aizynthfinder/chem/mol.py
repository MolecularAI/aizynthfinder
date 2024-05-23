""" Module containing classes to deal with Molecules - mostly wrappers around rdkit routines.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

from aizynthfinder.utils.bonds import sort_bonds
from aizynthfinder.utils.exceptions import MoleculeException

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Callable,
        Dict,
        List,
        Optional,
        RdMol,
        Sequence,
        Tuple,
        Union,
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
        self,
        rd_mol: Optional[RdMol] = None,
        smiles: Optional[str] = None,
        sanitize: bool = False,
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

    def fingerprint(
        self, radius: int, nbits: int = 2048, chiral: bool = False
    ) -> np.ndarray:
        """
        Returns the Morgan fingerprint of the molecule

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint
        :param chiral: if True, include chirality information
        :return: the fingerprint
        """
        key = radius, nbits

        if key not in self._fingerprints:
            self.sanitize()
            bitvect = AllChem.GetMorganFingerprintAsBitVect(
                self.rd_mol, *key, useChirality=chiral
            )
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

    def remove_atom_mapping(self, exceptions: Optional[Sequence[int]] = None) -> None:
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
        transform: Optional[int] = None,
        rd_mol: Optional[RdMol] = None,
        smiles: Optional[str] = None,
        sanitize: bool = False,
        mapping_update_callback: Optional[Callable[["TreeMolecule"], None]] = None,
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)
        self.parent = parent
        if transform is None and parent and parent.transform is not None:
            self.transform: int = parent.transform + 1
        else:
            self.transform = transform or 0

        self.original_smiles = smiles
        self.mapped_mol = Chem.Mol(self.rd_mol)
        self._atom_bonds: List[Tuple[int, int]] = []
        if not self.parent:
            self._set_atom_mappings()
        elif mapping_update_callback is not None:
            mapping_update_callback(self)

        AllChem.SanitizeMol(self.mapped_mol)
        self.mapped_smiles = Chem.MolToSmiles(self.mapped_mol)

        if self.parent:
            self.remove_atom_mapping()

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

    @property
    def mapped_atom_bonds(self) -> List[Tuple[int, int]]:
        """Return a list of atom bonds as tuples on the mapped atom indices"""
        bonds = []
        for bond in self.mapped_mol.GetBonds():
            bonds.append((bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()))

        self._atom_bonds = [
            (self.index_to_mapping[atom_index1], self.index_to_mapping[atom_index2])
            for atom_index1, atom_index2 in bonds
        ]
        return self._atom_bonds

    def get_bonds_in_molecule(
        self, query_bonds: Sequence[Sequence[int]]
    ) -> Sequence[Sequence[int]]:
        """
        Get bonds (from a list of bonds) that are present in the molecule.
        :param bonds: List of bond (atom pairs)
        :return: A list of bonds
        """
        molecule_bonds = sort_bonds(self.mapped_atom_bonds)
        query_bonds = sort_bonds(query_bonds)
        bonds_in_mol = [bond for bond in query_bonds if bond in molecule_bonds]
        return bonds_in_mol

    def has_all_focussed_bonds(self, bonds: Sequence[Sequence[int]]) -> bool:
        """Checks that the focussed bonds exist in the target molecule's atom bonds.

        :param bonds: Focussed bonds.
        :param target_mol: The target molecule.

        :return: A boolean indicating if the input bonds exist in the target molecule.
        """
        bonds_in_mol = self.get_bonds_in_molecule(bonds)
        return len(bonds_in_mol) == len(bonds)

    def _set_atom_mappings(self) -> None:
        atom_mappings = [
            atom.GetAtomMapNum()
            for atom in self.mapped_mol.GetAtoms()
            if atom.GetAtomMapNum() != 0
        ]

        mapper = max(atom_mappings) + 1 if atom_mappings else 1
        self._atom_mappings = {}
        for atom_index, atom in enumerate(self.mapped_mol.GetAtoms()):
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(mapper)
                mapper += 1
            self._atom_mappings[atom.GetAtomMapNum()] = atom_index


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
        self,
        rd_mol: Optional[RdMol] = None,
        smiles: Optional[str] = None,
        sanitize: bool = False,
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, _) -> bool:
        return False


def none_molecule() -> UniqueMolecule:
    """Return an empty molecule"""
    return UniqueMolecule(rd_mol=Chem.MolFromSmiles(""))
