"""
Tests to check the functionality of the bonds script.
"""
from aizynthfinder.chem import TreeMolecule, SmilesBasedRetroReaction
from aizynthfinder.utils.bonds import BrokenBonds


def test_focussed_bonds_broken():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    reaction = SmilesBasedRetroReaction(
        mol,
        mapped_prod_smiles="[CH3:1][NH:2][C:3](C)=[O:4]",
        reactants_str="C[C:3](=[O:4])O.[CH3:1][NH:2]",
    )
    focussed_bonds = [(1, 2), (3, 4), (3, 2)]
    broken_bonds = BrokenBonds(focussed_bonds)
    broken_focussed_bonds = broken_bonds(reaction)

    assert broken_focussed_bonds == [(2, 3)]


def test_focussed_bonds_not_broken():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    reaction = SmilesBasedRetroReaction(
        mol,
        mapped_prod_smiles="[CH3:1][NH:2][C:3](C)=[O:4]",
        reactants_str="C[C:3](=[O:4])O.[CH3:1][NH:2]",
    )
    focussed_bonds = [(1, 2), (3, 4)]
    broken_bonds = BrokenBonds(focussed_bonds)
    broken_focussed_bonds = broken_bonds(reaction)

    assert broken_focussed_bonds == []
    assert mol.has_all_focussed_bonds(focussed_bonds) is True


def test_focussed_bonds_not_in_target_mol():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    focussed_bonds = [(1, 4)]

    assert mol.has_all_focussed_bonds(focussed_bonds) is False
