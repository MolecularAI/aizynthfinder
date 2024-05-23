""" Module containing classes to deal with Reactions.
"""
from __future__ import annotations

import abc
import hashlib
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
from rdchiral import main as rdc

try:
    from rdchiral.bonds import get_atoms_across_double_bonds
    from rdchiral.initialization import BondDirOpposite
except ImportError:
    RDCHIRAL_CPP = True
else:
    RDCHIRAL_CPP = False
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondDir, BondStereo, ChiralType

from aizynthfinder.chem.mol import (
    Molecule,
    MoleculeException,
    TreeMolecule,
)
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.chem.mol import UniqueMolecule
    from aizynthfinder.utils.type_utils import (
        Any,
        Iterable,
        List,
        Optional,
        RdReaction,
        Set,
        StrDict,
        Tuple,
        Union,
    )

if RDCHIRAL_CPP:
    logger().warning(
        "WARNING: C++ version of RDChiral is supported, but with limited functionality"
    )


class _ReactionInterfaceMixin:
    """
    Mixin class to define a common interface for all reaction class

    The methods `_products_getter` and `_reactants_getter` needs to be implemented by subclasses
    """

    def fingerprint(
        self, radius: int, nbits: Optional[int] = None, chiral: bool = False
    ) -> np.ndarray:
        """
        Returns a difference fingerprint

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint. If not given it will use RDKit default, defaults to None
        :param chiral: if True, include chirality information
        :return: the fingerprint
        """
        product_fp = sum(
            mol.fingerprint(radius, nbits, chiral) for mol in self._products_getter()  # type: ignore
        )
        reactants_fp = sum(
            mol.fingerprint(radius, nbits, chiral) for mol in self._reactants_getter()  # type: ignore
        )
        return reactants_fp - product_fp  # type: ignore

    def hash_list(self) -> List[str]:
        """
        Return all the products and reactants as hashed SMILES

        :return: the hashes of the SMILES string
        """
        mols = self.reaction_smiles().replace(".", ">>").split(">>")
        return [hashlib.sha224(mol.encode("utf8")).hexdigest() for mol in mols]

    def hash_key(self) -> str:
        """
        Return a code that can be use to identify the reaction

        :return: the hash code
        """
        reactants = sorted([mol.inchi_key for mol in self._reactants_getter()])  # type: ignore
        products = sorted([mol.inchi_key for mol in self._products_getter()])  # type: ignore
        hash_ = hashlib.sha224()
        for item in reactants + [">>"] + products:
            hash_.update(item.encode())
        return hash_.hexdigest()

    def rd_reaction_from_smiles(self) -> RdReaction:
        """
        The reaction as a RDkit reaction object but created from the reaction smiles
        instead of the SMARTS of the template.

        :return: the reaction object
        """
        return AllChem.ReactionFromSmarts(self.reaction_smiles(), useSmiles=True)

    def reaction_smiles(self) -> str:
        """
        Get the reaction SMILES, i.e. the SMILES of the reactants and products joined together

        :return: the SMILES
        """
        reactants = ".".join(mol.smiles for mol in self._reactants_getter())  # type: ignore
        products = ".".join(mol.smiles for mol in self._products_getter())  # type: ignore
        return f"{reactants}>>{products}"


class RetroReaction(abc.ABC, _ReactionInterfaceMixin):
    """
    A retrosynthesis reaction. Only a single molecule is the reactant.

    This is an abstract class and child classes needs to implement the `_apply` and `_make_smiles` functions
    that should create the reactants molecule objects and the reaction SMILES representation, respectively.

    :ivar mol: the TreeMolecule object that this reaction is applied to
    :ivar index: a unique index of this reaction,
                 to count for the fact that a reaction can produce more than one outcome
    :ivar metadata: meta data associated with the reaction

    :param mol: the molecule
    :param index: the index, defaults to 0
    :param metadata: some meta data
    :params kwargs: any extra parameters for child classes
    """

    _required_kwargs: List[str] = []

    def __init__(
        self,
        mol: TreeMolecule,
        index: int = 0,
        metadata: Optional[StrDict] = None,
        **kwargs: Any,
    ) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise KeyError(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self.mol = mol
        self.index = index
        self.metadata: StrDict = metadata or {}
        self._reactants: Optional[Tuple[Tuple[TreeMolecule, ...], ...]] = None
        self._smiles: Optional[str] = None
        self._kwargs: StrDict = kwargs

    @classmethod
    def from_serialization(
        cls, init_args: StrDict, reactants: List[List[TreeMolecule]]
    ) -> "RetroReaction":
        """
        Create an object from a serialization. It does
        1) instantiate an object using the `init_args` and
        2) set the reactants to a tuple-form of `reactants

        :param init_args: the arguments passed to the `__init__` method
        :param reactants: the reactants
        :return: the deserialized object
        """
        obj = cls(**init_args)
        obj._reactants = tuple(tuple(mol for mol in lst_) for lst_ in reactants)
        return obj

    def __str__(self) -> str:
        return f"reaction on molecule {self.mol.smiles}"

    @property
    def reactants(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        """
        Returns the reactant molecules.
        Apply the reaction if necessary.

        :return: the products of the reaction
        """
        if not self._reactants:
            self._reactants = self._apply()
        return self._reactants

    @property
    def smiles(self) -> str:
        """
        The reaction as a SMILES

        :return: the SMILES
        """
        if self._smiles is None:
            try:
                self._smiles = self._make_smiles()
            except ValueError:
                self._smiles = ""  # noqa
        return self._smiles

    @property
    def unqueried(self) -> bool:
        """
        Return True if the reactants has never been retrieved
        """
        return self._reactants is None

    def copy(self, index: Optional[int] = None) -> "RetroReaction":
        """
        Shallow copy of this instance.

        :param index: new index, defaults to None
        :return: the copy
        """
        # pylint: disable=protected-access
        index = index if index is not None else self.index
        new_reaction = self.__class__(
            self.mol, index, dict(self.metadata), **self._kwargs
        )
        new_reaction._reactants = tuple(mol_list for mol_list in self._reactants or [])
        new_reaction._smiles = self._smiles
        return new_reaction

    def mapped_reaction_smiles(self) -> str:
        """
        Get the mapped reaction SMILES if it exists
        :return: the SMILES
        """
        reactants = self.mol.mapped_smiles
        products = ".".join(mol.mapped_smiles for mol in self._products_getter())
        return reactants + ">>" + products

    def to_dict(self) -> StrDict:
        """
        Return the retro reaction as dictionary
        This dictionary is not suitable for serialization, but is used by other serialization routines
        The elements of the dictionary can be used to instantiate a new reaction object
        """
        return {
            "mol": self.mol,
            "index": self.index,
            "metadata": dict(self.metadata),
        }

    @abc.abstractmethod
    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        pass

    @abc.abstractmethod
    def _make_smiles(self) -> str:
        pass

    def _products_getter(self) -> Tuple[TreeMolecule, ...]:
        return self.reactants[self.index]

    def _reactants_getter(self) -> List[TreeMolecule]:
        return [self.mol]

    @staticmethod
    def _update_unmapped_atom_num(mol: TreeMolecule, exclude_nums: Set[int]) -> None:
        mapped_nums = {num for num in mol.mapping_to_index.keys() if 0 < num < 900}
        offset = max(mapped_nums) + 1 if mapped_nums else 1
        for atom in mol.mapped_mol.GetAtoms():
            if 0 < atom.GetAtomMapNum() < 900:
                continue
            while offset in exclude_nums:
                offset += 1
            atom.SetAtomMapNum(offset)
            exclude_nums.add(offset)


class TemplatedRetroReaction(RetroReaction):
    """
    A retrosynthesis reaction that uses a reaction SMARTS and RDChiral to produce reactant molecules.
    The SMILES representation of the reaction is the SMARTS (modified by RDKit)

    :param mol: the molecule
    :param index: the index, defaults to 0
    :param metadata: some meta data
    :param smarts: a string representing the template
    """

    _required_kwargs = ["smarts"]

    def __init__(
        self,
        mol: TreeMolecule,
        index: int = 0,
        metadata: Optional[StrDict] = None,
        **kwargs: Any,
    ):
        super().__init__(mol, index, metadata, **kwargs)
        self.smarts: str = kwargs["smarts"]
        self._use_rdchiral: bool = kwargs.get("use_rdchiral", True)
        self._rd_reaction: Optional[RdReaction] = None

    def __str__(self) -> str:
        return (
            f"retro reaction from template {self.smarts} on molecule {self.mol.smiles}"
        )

    @property
    def rd_reaction(self) -> RdReaction:
        """Return the RDKit reaction created from the SMART"""
        if self._rd_reaction is None:
            self._rd_reaction = AllChem.ReactionFromSmarts(self.smarts)
        return self._rd_reaction

    def to_dict(self) -> StrDict:
        dict_ = super().to_dict()
        dict_["smarts"] = self.smarts
        return dict_

    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        if self._use_rdchiral:
            return self._apply_with_rdchiral()
        return self._apply_with_rdkit()

    def _apply_with_rdchiral(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        """
        Apply a reactions smarts to a molecule and return the products (reactants for retro templates)
        Will try to sanitize the reactants, and if that fails it will not return that molecule
        """
        reaction = rdc.rdchiralReaction(self.smarts)
        rct = _RdChiralProductWrapper(self.mol)
        try:
            reactants = rdc.rdchiralRun(reaction, rct, keep_mapnums=True)
        except RuntimeError as err:
            logger().debug(
                f"Runtime error in RDChiral with template {self.smarts} on {self.mol.smiles}\n{err}"
            )
            reactants = []
        except KeyError as err:
            logger().debug(
                f"Index error in RDChiral with template {self.smarts} on {self.mol.mapped_smiles}\n{err}"
            )
            reactants = []

        # Turning rdchiral outcome into rdkit tuple of tuples to maintain compatibility
        outcomes = []
        for reactant_str in reactants:
            smiles_list = reactant_str.split(".")
            exclude_nums = set(self.mol.mapping_to_index.keys())
            update_func = partial(
                self._update_unmapped_atom_num, exclude_nums=exclude_nums
            )
            try:
                rct_objs = tuple(
                    TreeMolecule(
                        parent=self.mol,
                        smiles=smi,
                        sanitize=True,
                        mapping_update_callback=update_func,
                    )
                    for smi in smiles_list
                )
            except MoleculeException:
                pass
            else:
                outcomes.append(rct_objs)
        self._reactants = tuple(outcomes)

        return self._reactants

    def _apply_with_rdkit(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        rxn = AllChem.ReactionFromSmarts(self.smarts)
        try:
            reactants_list = rxn.RunReactants([self.mol.mapped_mol])
        except:  # pylint: disable=bare-except
            reactants_list = []

        outcomes = []
        for reactants in reactants_list:
            exclude_nums = set(self.mol.mapping_to_index.keys())
            update_func = partial(self._inherit_atom_mapping, exclude_nums=exclude_nums)
            try:
                mols = tuple(
                    TreeMolecule(
                        parent=self.mol,
                        rd_mol=mol,
                        sanitize=True,
                        mapping_update_callback=update_func,
                    )
                    for mol in reactants
                )
            except MoleculeException:
                pass
            else:
                outcomes.append(mols)
        self._reactants = tuple(outcomes)

        return self._reactants

    def _make_smiles(self):
        return AllChem.ReactionToSmiles(self.rd_reaction)

    def _inherit_atom_mapping(self, mol: TreeMolecule, exclude_nums: Set[int]) -> None:
        """
        Update the internal atom mapping dictionary by inspecting the `reaction_atom_idx`
        property of the atoms and comparing it with the parent-molecule.

        This is used for child molecules created by RDKit reaction application.
        RDChiral takes care of this automatically.
        """
        if mol.parent is None:
            return

        for atom in mol.mapped_mol.GetAtoms():
            if not atom.HasProp("react_atom_idx"):
                continue
            index = atom.GetProp("react_atom_idx")
            mapping = mol.parent.index_to_mapping.get(int(index))
            if mapping:
                atom.SetAtomMapNum(mapping)

        self._update_unmapped_atom_num(mol, exclude_nums)


class SmilesBasedRetroReaction(RetroReaction):
    """
    A retrosynthesis reaction where the SMILES of the reactants are given on initiation

    The SMILES representation of the reaction is the reaction SMILES

    :param mol: the molecule
    :param index: the index, defaults to 0
    :param metadata: some meta data
    :param reactants_str: a dot-separated string of reactant SMILES strings
    """

    _required_kwargs = ["reactants_str"]

    def __init__(
        self,
        mol: TreeMolecule,
        index: int = 0,
        metadata: Optional[StrDict] = None,
        **kwargs: Any,
    ):
        super().__init__(mol, index, metadata, **kwargs)
        self.reactants_str: str = kwargs["reactants_str"]
        self._mapped_prod_smiles = kwargs.get("mapped_prod_smiles")

    def __str__(self) -> str:
        return (
            f"retro reaction on molecule {self.mol.smiles} giving {self.reactants_str}"
        )

    def to_dict(self) -> StrDict:
        dict_ = super().to_dict()
        dict_["reactants_str"] = self.reactants_str
        dict_["mapped_prod_smiles"] = self._mapped_prod_smiles
        return dict_

    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        outcomes = []
        smiles_list = self.reactants_str.split(".")

        exclude_nums = set(self.mol.mapping_to_index.keys())
        update_func = partial(self._remap, exclude_nums=exclude_nums)
        try:
            rct = tuple(
                TreeMolecule(
                    parent=self.mol,
                    smiles=smi,
                    sanitize=True,
                    mapping_update_callback=update_func,
                )
                for smi in smiles_list
            )
        except MoleculeException:
            pass
        else:
            outcomes.append(rct)
        self._reactants = tuple(outcomes)

        return self._reactants

    def _remap(self, mol: TreeMolecule, exclude_nums: Set[int]) -> None:
        """Find the mapping between parent and child and then re-map the child molecule"""
        if not self._mapped_prod_smiles:
            self._update_unmapped_atom_num(mol, exclude_nums)
            return

        parent_remapping = {}
        pmol = Molecule(smiles=self._mapped_prod_smiles, sanitize=True)
        for atom_idx1, atom_idx2 in enumerate(
            pmol.rd_mol.GetSubstructMatch(self.mol.mapped_mol)
        ):
            atom1 = self.mol.mapped_mol.GetAtomWithIdx(atom_idx1)
            atom2 = pmol.rd_mol.GetAtomWithIdx(atom_idx2)
            if atom1.GetAtomMapNum() > 0 and atom2.GetAtomMapNum() > 0:
                parent_remapping[atom2.GetAtomMapNum()] = atom1.GetAtomMapNum()

        for atom in mol.mapped_mol.GetAtoms():
            if atom.GetAtomMapNum() and atom.GetAtomMapNum() in parent_remapping:
                atom.SetAtomMapNum(parent_remapping[atom.GetAtomMapNum()])
            else:
                atom.SetAtomMapNum(0)

        self._update_unmapped_atom_num(mol, exclude_nums)

    def _make_smiles(self):
        rstr = ".".join(reactant.smiles for reactant in self.reactants[0])
        return f"{self.mol.smiles}>>{rstr}"


class FixedRetroReaction(_ReactionInterfaceMixin):
    """
    A retrosynthesis reaction that has the same interface as `RetroReaction`
    but it is fixed so it does not support SMARTS application or any creation of reactants.

    The reactants are set by using the `reactants` property.

    :ivar mol: the UniqueMolecule object that this reaction is applied to
    :ivar smiles: the SMILES representation of the RDKit reaction
    :ivar metadata: meta data associated with the reaction
    :ivar reactants: the reactants of this reaction

    :param mol: the molecule
    :param smiles: the SMILES of the reaction
    :param metadata: some meta data
    """

    def __init__(
        self,
        mol: UniqueMolecule,
        smiles: str = "",
        metadata: Optional[StrDict] = None,
    ) -> None:
        self.mol = mol
        self.smiles = smiles
        self.metadata = metadata or {}
        self.reactants: Tuple[Tuple[UniqueMolecule, ...], ...] = ()

    def copy(self) -> "FixedRetroReaction":
        """
        Shallow copy of this instance.

        :return: the copy
        """
        new_reaction = FixedRetroReaction(self.mol, self.smiles, self.metadata)
        new_reaction.reactants = tuple(mol_list for mol_list in self.reactants)
        return new_reaction

    def to_smiles_based_retroreaction(self) -> SmilesBasedRetroReaction:
        """
        Convert a FixedRetroReaction to a SmilesBasedRetroReaction.

        :return: the SmilesBasedRetroReaction.
        """
        if self.metadata and "mapped_reaction_smiles" in self.metadata.keys():
            mapped_reaction_smiles = self.metadata["mapped_reaction_smiles"]
        else:
            mapped_reaction_smiles = self.reaction_smiles()

        mapped_reaction_smiles = mapped_reaction_smiles.split(">>")
        product = mapped_reaction_smiles[0]
        reactants = mapped_reaction_smiles[1]

        return SmilesBasedRetroReaction(
            mol=TreeMolecule(smiles=product, parent=None),
            mapped_prod_smiles=product,
            reactants_str=reactants,
        )

    def _products_getter(self) -> Tuple[UniqueMolecule, ...]:
        return self.reactants[0]

    def _reactants_getter(self) -> List[UniqueMolecule]:
        return [self.mol]


def hash_reactions(
    reactions: Union[Iterable[RetroReaction], Iterable[FixedRetroReaction]],
    sort: bool = True,
) -> str:
    """
    Creates a hash for a list of reactions

    :param reactions: the reactions to hash
    :param sort: if True will sort all molecules, defaults to True
    :return: the hash string
    """
    hash_list = []
    for reaction in reactions:
        hash_list.extend(reaction.hash_list())
    if sort:
        hash_list.sort()
    hash_list_str = ".".join(hash_list)
    return hashlib.sha224(hash_list_str.encode("utf8")).hexdigest()


class _RdChiralProductWrapper:
    """
    Reimplementation of `rdchiralReaction`
    to preserve product molecule already created
    """

    # pylint: disable=W0106,C0103
    def __init__(self, product: TreeMolecule) -> None:
        product.sanitize()
        self.reactant_smiles = product.smiles

        # Initialize into RDKit mol
        self.reactants = Chem.Mol(product.mapped_mol.ToBinary())
        Chem.AssignStereochemistry(self.reactants, flagPossibleStereoCenters=True)
        self.reactants.UpdatePropertyCache(strict=False)

        self.atoms_r = {a.GetAtomMapNum(): a for a in self.reactants.GetAtoms()}
        self.idx_to_mapnum = lambda idx: self.reactants.GetAtomWithIdx(
            idx
        ).GetAtomMapNum()

        # Create copy of molecule without chiral information, used with
        # RDKit's naive runReactants
        self.reactants_achiral = Chem.Mol(product.rd_mol.ToBinary())
        [
            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
            for a in self.reactants_achiral.GetAtoms()
        ]
        [
            (b.SetStereo(BondStereo.STEREONONE), b.SetBondDir(BondDir.NONE))
            for b in self.reactants_achiral.GetBonds()
        ]

        # Pre-list reactant bonds (for stitching broken products)
        self.bonds_by_mapnum = [
            (b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum(), b)
            for b in self.reactants.GetBonds()
        ]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.bond_dirs_by_mapnum = {}
        for i, j, b in self.bonds_by_mapnum:
            if b.GetBondDir() != BondDir.NONE:
                self.bond_dirs_by_mapnum[(i, j)] = b.GetBondDir()
                self.bond_dirs_by_mapnum[(j, i)] = BondDirOpposite[b.GetBondDir()]

        # Get atoms across double bonds defined by mapnum
        self.atoms_across_double_bonds = get_atoms_across_double_bonds(self.reactants)


if RDCHIRAL_CPP:

    def _wrapper(mol):
        return rdc.rdchiralReactants(mol.mapped_smiles)

    _RdChiralProductWrapper = _wrapper  # type: ignore
