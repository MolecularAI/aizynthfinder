import pickle

import pytest
from rdkit import Chem

from aizynthfinder.utils.sc_score import SCScore

# Dummy, tiny SCScore model
_weights0 = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]
_weights1 = [[1, 1], [1, 1]]
_weights = [_weights0, _weights1]
_biases = [[0, 0], [0]]


def test_scscore(tmpdir):
    filename = str(tmpdir / "dummy.pickle")
    with open(filename, "wb") as fileobj:
        pickle.dump((_weights, _biases), fileobj)
    scorer = SCScore(filename, 5)
    mol = Chem.MolFromSmiles("C")

    assert pytest.approx(scorer(mol), abs=1e-3) == 4.523
