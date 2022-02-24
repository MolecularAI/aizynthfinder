import numpy as np

import pytest

from aizynthfinder.context.cost import MoleculeCost
from aizynthfinder.search.retrostar.cost import RetroStarCost
from aizynthfinder.chem import Molecule


def test_retrostar_cost(setup_mocked_model):
    mol = Molecule(smiles="CCCC")

    cost = RetroStarCost(model_path="dummy", fingerprint_length=10, dropout_rate=0.0)
    assert pytest.approx(cost(mol), abs=0.001) == 30


def test_load_cost_from_config(setup_mocked_model):
    cost = MoleculeCost()

    dict_ = {
        "aizynthfinder.search.retrostar.cost.RetroStarCost": {"model_path": "dummy"}
    }
    cost.load_from_config(**dict_)

    assert len(cost) == 2
