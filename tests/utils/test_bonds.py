"""
Tests to check the functionality of the bonds script.
"""
from aizynthfinder.chem import TreeMolecule, SmilesBasedRetroReaction
from aizynthfinder.utils.bonds import BrokenBonds


def test_focused_bonds_broken():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    reaction = SmilesBasedRetroReaction(
        mol,
        mapped_prod_smiles="[CH3:1][NH:2][C:3](C)=[O:4]",
        reactants_str="C[C:3](=[O:4])O.[CH3:1][NH:2]",
    )
    focused_bonds = [(1, 2), (3, 4), (3, 2)]
    broken_bonds = BrokenBonds(focused_bonds)
    broken_focused_bonds = broken_bonds(reaction)

    assert broken_focused_bonds == [(2, 3)]


def test_focused_bonds_not_broken():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    reaction = SmilesBasedRetroReaction(
        mol,
        mapped_prod_smiles="[CH3:1][NH:2][C:3](C)=[O:4]",
        reactants_str="C[C:3](=[O:4])O.[CH3:1][NH:2]",
    )
    focused_bonds = [(1, 2), (3, 4)]
    broken_bonds = BrokenBonds(focused_bonds)
    broken_focused_bonds = broken_bonds(reaction)

    assert broken_focused_bonds == []
    assert mol.has_all_focused_bonds(focused_bonds) is True


def test_focused_bonds_not_in_target_mol():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    focused_bonds = [(1, 4)]

    assert mol.has_all_focused_bonds(focused_bonds) is False
