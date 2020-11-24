import pytest
from rdkit import Chem

from aizynthfinder.chem import MoleculeException, Molecule, TreeMolecule


def test_no_input():
    with pytest.raises(MoleculeException):
        Molecule()


def test_create_with_mol():
    rd_mol = Chem.MolFromSmiles("O")

    mol = Molecule(rd_mol=rd_mol)

    assert mol.smiles == "O"


def test_create_with_smiles():
    mol = Molecule(smiles="O")

    assert Chem.MolToSmiles(mol.rd_mol) == "O"


def test_inchi_key():
    mol = Molecule(smiles="O")

    assert mol.inchi_key == "XLYOFNOQVPJJNP-UHFFFAOYSA-N"


def test_fingerprint():
    mol = Molecule(smiles="O")

    assert sum(mol.fingerprint(2)) == 1

    assert sum(mol.fingerprint(2, 10)) == 1


def test_sanitize():
    mol = Molecule(smiles="O", sanitize=True)

    assert Chem.MolToSmiles(mol.rd_mol) == "O"

    mol = Molecule(smiles="c1ccccc1(C)(C)")

    with pytest.raises(MoleculeException):
        mol.sanitize()


def test_equality():
    mol1 = Molecule(smiles="CCCCO")
    mol2 = Molecule(smiles="OCCCC")

    assert mol1 == mol2


def test_reaction(simple_actions):
    mol = TreeMolecule(parent=None, smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
    reactions, _ = simple_actions(mol)

    products1 = reactions[0].apply()
    assert products1[0][0].smiles == "CCCCOc1ccc(CC(=O)Cl)cc1"
    assert products1[0][1].smiles == "CNO"

    products2 = reactions[2].apply()
    assert products2 == ()
