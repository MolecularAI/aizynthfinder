from aizynthfinder.context.cost import (
    ZeroCost,
    MoleculeCost,
)
from aizynthfinder.chem import Molecule


def test_zero_cost():
    assert ZeroCost()(None) == 0.0


def test_create_mol_cost():
    mol = Molecule(smiles="CCCC")

    assert MoleculeCost()(mol) == 0.0


def test_mol_cost_cache(mocker):
    mol = Molecule(smiles="CCCC")
    mocker.patch("aizynthfinder.context.cost.ZeroCost.__call__")
    cost = MoleculeCost()
    cost[cost.selection].__call__.return_value = 5.0

    assert cost(mol) == 5.0
    assert cost(mol) == 5.0
    cost[cost.selection].__call__.assert_called_once()


def test_load_cost():
    cost = MoleculeCost()

    cost.load(ZeroCost())
    assert len(cost) == 1


def test_load_cost_from_config():
    cost = MoleculeCost()

    cost.load_from_config(zero={})

    assert len(cost) == 1
