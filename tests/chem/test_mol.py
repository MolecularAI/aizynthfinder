import pytest
from rdkit import Chem

from aizynthfinder.chem import MoleculeException, Molecule


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


def test_inchi():
    mol = Molecule(smiles="O")

    assert mol.inchi == "InChI=1S/H2O/h1H2"


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

    mol.sanitize(raise_exception=False)
    assert mol.smiles == "CC1(C)CCCCC1"


def test_equality():
    mol1 = Molecule(smiles="CCCCO")
    mol2 = Molecule(smiles="OCCCC")

    assert mol1 == mol2


def test_basic_equality():
    mol1 = Molecule(smiles="CC[C@@H](C)O")  # R-2-butanol
    mol2 = Molecule(smiles="CC[C@H](C)O")  # S-2-butanol

    assert mol1 != mol2
    assert mol1.basic_compare(mol2)


def test_has_atom_mapping():
    mol1 = Molecule(smiles="CCCCO")
    mol2 = Molecule(smiles="C[C:5]CCO")

    assert not mol1.has_atom_mapping()
    assert mol2.has_atom_mapping()


def test_remove_atom_mapping():
    mol = Molecule(smiles="C[C:5]CCO")

    assert mol.has_atom_mapping()

    mol.remove_atom_mapping()

    assert not mol.has_atom_mapping()
