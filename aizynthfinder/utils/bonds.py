"""Module containing a class to identify broken focussed bonds
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import List, Sequence, Tuple
    from aizynthfinder.chem.mol import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction


class BrokenBonds:
    """
    A class to keep track of focussed bonds breaking in a target molecule.

    :param focussed_bonds: A list of focussed bond pairs. The bond pairs are represented
        as tuples of size 2. These bond pairs exist in the target molecule's atom bonds.
    """

    def __init__(self, focussed_bonds: Sequence[Sequence[int]]) -> None:
        self.focussed_bonds = sort_bonds(focussed_bonds)
        self.filtered_focussed_bonds: List[Tuple[int, int]] = []

    def __call__(self, reaction: RetroReaction) -> List[Tuple[int, int]]:
        """
        Provides a list of focussed bonds that break in any of the molecule's reactants.

        :param reaction: A retro reaction.
        :return: A list of all the focussed bonds that broke within the reactants
            that constitute the target molecule.
        """
        self.filtered_focussed_bonds = self._get_filtered_focussed_bonds(reaction.mol)
        if not self.filtered_focussed_bonds:
            return []

        molecule_bonds = []
        for reactant in reaction.reactants[reaction.index]:
            molecule_bonds += reactant.mapped_atom_bonds

        broken_focussed_bonds = self._get_broken_frozen_bonds(
            sort_bonds(molecule_bonds)
        )
        return broken_focussed_bonds

    def _get_broken_frozen_bonds(
        self,
        molecule_bonds: List[Tuple[int, int]],
    ) -> List[Tuple[int, int]]:
        broken_focussed_bonds = list(
            set(self.filtered_focussed_bonds) - set(molecule_bonds)
        )
        return broken_focussed_bonds

    def _get_filtered_focussed_bonds(
        self, molecule: TreeMolecule
    ) -> List[Tuple[int, int]]:
        molecule_bonds = molecule.mapped_atom_bonds
        atom_maps = [atom_map for bonds in molecule_bonds for atom_map in bonds]

        filtered_focussed_bonds = []
        for idx1, idx2 in self.focussed_bonds:
            if idx1 in atom_maps and idx2 in atom_maps:
                filtered_focussed_bonds.append((idx1, idx2))
        return filtered_focussed_bonds


def sort_bonds(bonds: Sequence[Sequence[int]]) -> List[Tuple[int, int]]:
    return [tuple(sorted(bond)) for bond in bonds]  # type: ignore
