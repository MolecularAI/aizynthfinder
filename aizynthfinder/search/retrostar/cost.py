""" Module containing Retro* cost model """
from __future__ import annotations
from typing import TYPE_CHECKING
import pickle

import numpy as np

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Tuple, List
    from aizynthfinder.chem import Molecule


class RetroStarCost:
    """
    Encapsulation of a the original Retro* molecular cost model

    Numpy implementation of original pytorch model

    The predictions of the score is made on a Molecule object

    .. code-block::

        mol = Molecule(smiles="CCC")
        scorer = RetroStarCost()
        score = scorer(mol)

    The model provided when creating the scorer object should be a pickled
    tuple.
    The first item of the tuple should be a list of the model weights for each layer.
    The second item of the tuple should be a list of the model biases for each layer.

    :param model_path: the filename of the model weights and biases
    :param fingerprint_length: the number of bits in the fingerprint
    :param fingerprint_radius: the radius of the fingerprint
    :param dropout_rate: the dropout_rate
    """

    def __init__(
        self,
        model_path: str,
        fingerprint_length: int = 2048,
        fingerprint_radius: int = 2,
        dropout_rate: float = 0.1,
    ):
        self._dropout_prob = 1.0 - dropout_rate
        self._fingerprint_length = fingerprint_length
        self._fingerprint_radius = fingerprint_radius
        self._weights, self._biases = self._load_model(model_path)

    def __call__(self, mol: Molecule) -> float:
        # pylint: disable=invalid-name
        mol.sanitize()
        vec = mol.fingerprint(
            radius=self._fingerprint_radius, nbits=self._fingerprint_length
        )
        for W, b in zip(self._weights[:-1], self._biases[:-1]):
            vec = np.matmul(vec, W) + b
            vec *= vec > 0  # ReLU
            # Drop-out
            vec *= np.random.binomial(1, self._dropout_prob, size=vec.shape) / (
                self._dropout_prob
            )
        vec = np.matmul(vec, self._weights[-1]) + self._biases[-1]
        return float(np.log(1 + np.exp(vec)))

    def __repr__(self) -> str:
        return "retrostar"

    @staticmethod
    def _load_model(model_path: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:

        with open(model_path, "rb") as fileobj:
            weights, biases = pickle.load(fileobj)

        return (
            [np.asarray(item) for item in weights],
            [np.asarray(item) for item in biases],
        )
