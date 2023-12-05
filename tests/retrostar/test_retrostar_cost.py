import pytest

from aizynthfinder.chem import Molecule
from aizynthfinder.search.retrostar.cost import (
    MoleculeCost,
    RetroStarCost,
    ZeroMoleculeCost,
)


def test_retrostar_cost(setup_mocked_model):
    mol = Molecule(smiles="CCCC")

    cost = RetroStarCost(model_path="dummy", fingerprint_length=10, dropout_rate=0.0)
    assert pytest.approx(cost.calculate(mol), abs=0.001) == 30


def test_zero_molecule_cost():
    mol = Molecule(smiles="CCCC")

    cost = ZeroMoleculeCost().calculate(mol)
    assert cost == 0.0


def test_molecule_cost_zero(default_config):
    default_config.search.algorithm_config["molecule_cost"] = {
        "cost": "ZeroMoleculeCost"
    }
    mol = Molecule(smiles="CCCC")

    molecule_cost = MoleculeCost(default_config)(mol)
    assert molecule_cost == 0.0


def test_molecule_cost_retrostar(default_config, setup_mocked_model):
    default_config.search.algorithm_config["molecule_cost"] = {
        "cost": "aizynthfinder.search.retrostar.cost.RetroStarCost",
        "model_path": "dummy",
        "fingerprint_length": 10,
        "dropout_rate": 0.0,
    }
    mol = Molecule(smiles="CCCC")

    molecule_cost = MoleculeCost(default_config)(mol)
    assert pytest.approx(molecule_cost, abs=0.001) == 30
