""" Module containing classes to deal with Reactions.
"""
from __future__ import annotations
import hashlib
import abc
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.rdchem import KekulizeException
from rdchiral import main as rdc

from aizynthfinder.utils.logging import logger
from aizynthfinder.chem.mol import MoleculeException, Molecule, TreeMolecule

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Optional,
        Union,
        Tuple,
        List,
        RdReaction,
        StrDict,
        Iterable,
        Any,
    )
    from aizynthfinder.chem.mol import UniqueMolecule


class _ReactionInterfaceMixin:
    """
    Mixin class to define a common interface for all reaction class

    The methods `_products_getter` and `_reactants_getter` needs to be implemented by subclasses
    """

    def fingerprint(self, radius: int, nbits: int = None) -> np.ndarray:
        """
        Returns a difference fingerprint

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint. If not given it will use RDKit default, defaults to None
        :return: the fingerprint
        """
        product_fp = sum(
            mol.fingerprint(radius, nbits) for mol in self._products_getter()  # type: ignore
        )
        reactants_fp = sum(
            mol.fingerprint(radius, nbits) for mol in self._reactants_getter()  # type: ignore
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
        return "%s>>%s" % (reactants, products)


class Reaction(_ReactionInterfaceMixin):
    """
    An abstraction of a reaction. Encapsulate an RDKit reaction object and
    functions that can be applied to such a reaction.

    :ivar mols: the Molecule objects that this reaction are applied to
    :ivar smarts: the SMARTS representation of the reaction
    :ivar index: a unique index of this reaction,
                 to count for the fact that a reaction can produce more than one outcome
    :ivar metadata: meta data associated with the reaction

    :param mols: the molecules
    :param smarts: the SMARTS fragment
    :param index: the index, defaults to 0
    :param metadata: some meta data
    """

    def __init__(
        self,
        mols: List[Molecule],
        smarts: str,
        index: int = 0,
        metadata: StrDict = None,
    ) -> None:
        self.mols = mols
        self.smarts = smarts
        self.index = index
        self.metadata: StrDict = metadata or {}
        self._products: Optional[Tuple[Tuple[Molecule, ...], ...]] = None
        self._rd_reaction: Optional[RdReaction] = None
        self._smiles: Optional[str] = None

    @property
    def products(self) -> Tuple[Tuple[Molecule, ...], ...]:
        """
        Returns the product molecules.
        Apply the reaction if necessary.

        :return: the products of the reaction
        """
        if not self._products:
            self.apply()
            assert self._products is not None
        return self._products

    @property
    def rd_reaction(self) -> RdReaction:
        """
        The reaction as a RDkit reaction object

        :return: the reaction object
        """
        if not self._rd_reaction:
            self._rd_reaction = AllChem.ReactionFromSmarts(self.smarts)
        return self._rd_reaction

    @property
    def smiles(self) -> str:
        """
        The reaction as a SMILES

        :return: the SMILES
        """
        if self._smiles is None:
            try:
                self._smiles = AllChem.ReactionToSmiles(self.rd_reaction)
            except ValueError:
                self._smiles = ""  # noqa
        return self._smiles

    def apply(self) -> Tuple[Tuple[Molecule, ...], ...]:
        """
        Apply a reactions smarts to list of reactant and return the products

        Will try to sanitize the reactants, and if that fails it will not return that molecule

        :return: the products of the reaction
        """
        num_rectantant_templates = self.rd_reaction.GetNumReactantTemplates()
        reactants = tuple(mol.rd_mol for mol in self.mols[:num_rectantant_templates])
        products_list = self.rd_reaction.RunReactants(reactants)

        outcomes = []
        for products in products_list:
            try:
                mols = tuple(Molecule(rd_mol=mol, sanitize=True) for mol in products)
            except MoleculeException:
                pass
            else:
                outcomes.append(mols)
        self._products = tuple(outcomes)

        return self._products

    def _products_getter(self) -> Tuple[Molecule, ...]:
        return self.products[self.index]

    def _reactants_getter(self) -> List[Molecule]:
        return self.mols


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
        self, mol: TreeMolecule, index: int = 0, metadata: StrDict = None, **kwargs: Any
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

    def copy(self, index: int = None) -> "RetroReaction":
        """
        Shallow copy of this instance.

        :param index: new index, defaults to None
        :return: the copy
        """
        # pylint: disable=protected-access
        index = index if index is not None else self.index
        new_reaction = self.__class__(self.mol, index, self.metadata, **self._kwargs)
        new_reaction._reactants = tuple(mol_list for mol_list in self._reactants or [])
        new_reaction._smiles = self._smiles
        return new_reaction

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
        self, mol: TreeMolecule, index: int = 0, metadata: StrDict = None, **kwargs: Any
    ):
        super().__init__(mol, index, metadata, **kwargs)
        self.smarts: str = kwargs["smarts"]
        self._use_rdchiral: bool = kwargs.get("use_rdchiral", True)
        self._reaction_source: str = kwargs.get("reaction_source", "rdkit")
        self._template_fallback: str = kwargs.get("template_fallback", True)
        self.templates = kwargs.get("templates", None)
        self._rd_reaction: Optional[RdReaction] = None
        self.inchi_fail = kwargs.get("inchi_fail", 0)
        self.rdkit_fail = kwargs.get("rdkit_fail", 0)

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

    def forward_reaction(self) -> Reaction:
        """
        Create the forward reaction corresponding to the SMARTS and reactants of this reaction

        :return: the forward reaction
        """
        fwd_smarts = ">>".join(self.smarts.split(">>")[::-1])
        mols = [Molecule(rd_mol=mol.rd_mol) for mol in self.reactants[self.index]]
        return Reaction(mols=mols, smarts=fwd_smarts)

    def to_dict(self) -> StrDict:
        dict_ = super().to_dict()
        dict_["smarts"] = self.smarts
        return dict_

    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        if self._reaction_source == "rdkit":
            return self._apply_with_rdkit()
        elif self._reaction_source == "rdchiral":
            return self._apply_with_rdchiral()
        elif self._reaction_source == "templates":
            return self._apply_with_templates()

    def _apply_with_rdchiral(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        """
        Apply a reactions smarts to a molecule and return the products (reactants for retro templates)
        Will try to sanitize the reactants, and if that fails it will not return that molecule
        """
        reaction = rdc.rdchiralReaction(self.smarts)
        rct = rdc.rdchiralReactants(self.mol.smiles)
        try:
            reactants = rdc.rdchiralRun(reaction, rct)
        except RuntimeError as err:
            logger().debug(
                f"Runtime error in RDChiral with template {self.smarts} on {self.mol.smiles}\n{err}"
            )
            reactants = []
        except MoleculeException:
            if self._template_fallback: # TODO: come up with a better name
                return self._apply_with_templates()
            reactants = []

        # Turning rdchiral outcome into rdkit tuple of tuples to maintain compatibility
        outcomes = []
        for reactant_str in reactants:
            smiles_list = reactant_str.split(".")
            try:
                rct = tuple(
                    TreeMolecule(parent=self.mol, smiles=smi, sanitize=True)
                    for smi in smiles_list
                )
            except MoleculeException:
                pass
            else:
                outcomes.append(rct)
        self._reactants = tuple(outcomes)

        return self._reactants

    def _apply_with_rdkit(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        rxn = AllChem.ReactionFromSmarts(self.smarts)
        try:
            self.mol.sanitize()
        except MoleculeException:
            reactants_list = []
        try:
            reactants_list = rxn.RunReactants([self.mol.rd_mol])
        except:
            self.rdkit_fail += 1
            if self._template_fallback: # TODO: come up with a better name
                return self._apply_with_templates()
            reactants_list = []

        outcomes = []
        for reactants in reactants_list:
            try:
                mols = tuple(
                    TreeMolecule(parent=self.mol, rd_mol=mol, sanitize=True)
                    for mol in reactants
                )
            except MoleculeException:
                pass
            else:
                outcomes.append(mols)
        self._reactants = tuple(outcomes)

        return self._reactants
    
    def _apply_with_templates(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        # TODO: get column names from config file
        rxn = self.templates.loc[self.templates["retro_template"]==self.smarts]["retro_template"]
        reactants_list = rxn.values[0].split(">>")[0].split(".")
        rxn_fp = Chem.RDKFingerprint(self.mol.rd_mol)
        react_rp = [
            {"smarts": r,
             "mol": Chem.MolFromSmarts(r),
             "fp": Chem.RDKFingerprint(Chem.MolFromSmarts(r))}
            for r in reactants_list
            ]
        sims = [
            {"mol": x["mol"],
             "sim": FingerprintSimilarity(rxn_fp, x["fp"])}
            for x in react_rp
            ]
        sorted_sims = sorted(sims, key=lambda m: m["sim"], reverse=True)
        if len(sorted_sims) > 1:
            sorted_sims.pop(0)
        
        reactants = [r["mol"] for r in sorted_sims]
        reactants.append(self.mol.rd_mol)
        
        outcomes = []
        for mol in reactants:
            try:
                mols = (
                    TreeMolecule(parent=self.mol, rd_mol=mol,
                                 sanitize=True),
                )
                # Check for Inchi Conversion
                Chem.MolToInchi(mol)
            except (MoleculeException, KekulizeException):
                self.inchi_fail += 1
                pass
            else:
                outcomes.append(mols)

        self._reactants = tuple(outcomes)
        
        return self._reactants
    
    
    # def _apply_with_templates(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
    #     # TODO: get column names from config file
    #     rxn = self.templates.loc[self.templates["retro_template"]==self.smarts]
    #     reactants_list = []
    #     for r in rxn["reactants"].values:
    #         reactants_list.append((AllChem.MolFromSmiles(r),))
        
    #     outcomes = []
    #     for reactants in reactants_list:
    #         try:
    #             mols = tuple(
    #             TreeMolecule(parent=self.mol, rd_mol=mol, sanitize=True)
    #                 for mol in reactants
    #             )
    #         except MoleculeException:
    #             pass
    #         else:
    #             outcomes.append(mols)
    #     self._reactants = tuple(outcomes)
        
    #     return self._reactants        

    def _make_smiles(self):
        return AllChem.ReactionToSmiles(self.rd_reaction)


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
        self, mol: TreeMolecule, index: int = 0, metadata: StrDict = None, **kwargs: Any
    ):
        super().__init__(mol, index, metadata, **kwargs)
        self.reactants_str: str = kwargs["reactants_str"]

    def __str__(self) -> str:
        return (
            f"retro reaction on molecule {self.mol.smiles} giving {self.reactants_str}"
        )

    def to_dict(self) -> StrDict:
        dict_ = super().to_dict()
        dict_["reactants_str"] = self.reactants_str
        return dict_

    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        outcomes = []
        smiles_list = self.reactants_str.split(".")
        try:
            rct = tuple(
                TreeMolecule(parent=self.mol, smiles=smi, sanitize=True)
                for smi in smiles_list
            )
        except MoleculeException:
            pass
        else:
            outcomes.append(rct)
        self._reactants = tuple(outcomes)

        return self._reactants

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
        self, mol: UniqueMolecule, smiles: str = "", metadata: StrDict = None
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

    def _products_getter(self) -> Tuple[UniqueMolecule, ...]:
        return self.reactants[0]

    def _reactants_getter(self) -> List[UniqueMolecule]:
        return [self.mol]


def hash_reactions(
    reactions: Union[
        Iterable[Reaction], Iterable[RetroReaction], Iterable[FixedRetroReaction]
    ],
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
