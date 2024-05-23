""" Module contain a class that encapsulate the state of search tree node.
"""
from __future__ import annotations

import os
from typing import TYPE_CHECKING

from rdkit.Chem import Draw

from aizynthfinder.chem import TreeMolecule

if TYPE_CHECKING:
    from aizynthfinder.chem import MoleculeDeserializer, MoleculeSerializer
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Any, List, Optional, Sequence, StrDict

os.environ["PYTHONHASHSEED"] = "42"


class MctsState:
    """
    Encapsulation of the molecular state of a node.

    A state consists of an immutable list of molecules that are either solved
    (can be found in stock) or that potentially can be expanded to new molecules
    by applying a reaction on them.

    The class is hashable and comparable by the inchi keys of all the molecules.

    :ivar mols: the list of molecules
    :ivar expandable_mols: the list of molecules not in stock
    :ivar stock: the configured stock
    :ivar in_stock_list: for each molecule if they are in stock
    :ivar is_solved: is true if all molecules are in stock:
    :ivar max_transforms: the maximum of the transforms of the molecule
    :ivar is_terminal: is true if the all molecules are in stock or if the maximum transforms has been reached
    :ivar expandables_hash: an hash string computed on the expandable molecules

    :param mols: the molecules of the state
    :param config: settings of the tree search algorithm
    """

    def __init__(self, mols: Sequence[TreeMolecule], config: Configuration) -> None:
        self.mols = mols
        self.stock = config.stock
        self.in_stock_list = [mol in self.stock for mol in self.mols]
        self.expandable_mols = [
            mol for mol, in_stock in zip(self.mols, self.in_stock_list) if not in_stock
        ]
        self._stock_availability: Optional[List[str]] = None
        self.is_solved = all(self.in_stock_list)
        self.max_transforms = max(mol.transform for mol in self.mols)
        self.is_terminal = (
            self.max_transforms >= config.search.max_transforms
        ) or self.is_solved

        inchis = [mol.inchi_key for mol in self.mols]
        self._hash = hash(tuple(sorted(inchis)))

        inchis = [mol.inchi_key for mol in self.expandable_mols]
        self.expandables_hash = hash(tuple(sorted(inchis)))

    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MctsState):
            return False
        return self.__hash__() == other.__hash__()

    def __str__(self) -> str:
        """A string representation of the state (for print(state))"""
        string = (
            f"{str([mol.smiles for mol in self.mols])}\n"
            + f"{str([mol.transform for mol in self.mols])}\n"
            + f"{str([mol.parent.smiles if mol.parent else '-' for mol in self.mols])}\n"
            + f"{str(self.in_stock_list)}\n Solved: {self.is_solved}"
        )
        return string

    @classmethod
    def from_dict(
        cls, dict_: StrDict, config: Configuration, molecules: MoleculeDeserializer
    ) -> "MctsState":
        """
        Create a new state from a dictionary, i.e. deserialization

        :param dict_: the serialized state
        :type dict_: dict
        :param config: settings of the tree search algorithm
        :type config: Configuration
        :param molecules: the deserialized molecules
        :type molecules: MoleculeDeserializer
        :return: a deserialized state
        :rtype: State
        """
        mols = molecules.get_tree_molecules(dict_["mols"])
        return MctsState(mols, config)

    @property
    def stock_availability(self) -> List[str]:
        """
        Returns a list of availabilities for all molecules

        :return: the list
        :rtype: list of str
        """
        if not self._stock_availability:
            self._stock_availability = [
                self.stock.availability_string(mol) for mol in self.mols
            ]
        return self._stock_availability

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the state object to a dictionary

        :param molecule_store: the serialized molecules
        :type molecule_store: MolecularSerializer
        :return: the serialized state
        :rtype: dict
        """
        return {"mols": [molecule_store[mol] for mol in self.mols]}

    def to_image(self, ncolumns: int = 6) -> Any:
        """
        Constructs an image representation of the state

        :param ncolumns: number of molecules per row, defaults to 6
        :type ncolumns: int, optional
        :return: the image representation
        :rtype: a PIL image
        """
        for mol in self.mols:
            mol.sanitize()
        legends = self.stock_availability
        mols = [mol.rd_mol for mol in self.mols]
        return Draw.MolsToGridImage(mols, molsPerRow=ncolumns, legends=legends)
