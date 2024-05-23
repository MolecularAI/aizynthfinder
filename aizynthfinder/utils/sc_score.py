""" Module containing the implementation of the SC-score model for synthetic complexity scoring. """

from __future__ import annotations

import pickle
from typing import TYPE_CHECKING

import numpy as np
from rdkit.Chem import AllChem

from aizynthfinder.utils.math import (
    dense_layer_forward_pass,
    rectified_linear_unit,
    sigmoid,
)

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import RdMol, Sequence, Tuple


class SCScore:
    """
    Encapsulation of the SCScore model

    Re-write of the SCScorer from the scscorer package

    The predictions of the score is made with a sanitized instance of an RDKit molecule

    .. code-block::

        mol = Molecule(smiles="CCC", sanitize=True)
        scscorer = SCScorer("path_to_model")
        score = scscorer(mol.rd_mol)

    The model provided when creating the scorer object should be pickled tuple.
    The first item of the tuple should be a list of the model weights for each layer.
    The second item of the tuple should be a list of the model biases for each layer.

    :param model_path: the filename of the model weights and biases
    :param fingerprint_length: the number of bits in the fingerprint
    :param fingerprint_radius: the radius of the fingerprint
    """

    def __init__(
        self,
        model_path: str,
        fingerprint_length: int = 1024,
        fingerprint_radius: int = 2,
    ) -> None:
        self._fingerprint_length = fingerprint_length
        self._fingerprint_radius = fingerprint_radius
        self._weights, self._biases = self._load_model(model_path)
        self.score_scale = 5.0

    def __call__(self, rd_mol: RdMol) -> float:
        fingerprint = self._make_fingerprint(rd_mol)
        normalized_score = self.forward(fingerprint)
        sc_score = (1 + (self.score_scale - 1) * normalized_score)[0]
        return sc_score

    # pylint: disable=invalid-name
    def forward(self, x: np.ndarray) -> np.ndarray:
        """Forward pass with dense neural network"""
        for weights, bias in zip(self._weights[:-1], self._biases[:-1]):
            x = dense_layer_forward_pass(x, weights, bias, rectified_linear_unit)
        return dense_layer_forward_pass(x, self._weights[-1], self._biases[-1], sigmoid)

    def _load_model(
        self, model_path: str
    ) -> Tuple[Sequence[np.ndarray], Sequence[np.ndarray]]:
        """Returns neural network model parameters."""
        with open(model_path, "rb") as fileobj:
            weights, biases = pickle.load(fileobj)

        weights = [np.asarray(item) for item in weights]
        biases = [np.asarray(item) for item in biases]
        return weights, biases

    def _make_fingerprint(self, rd_mol: RdMol) -> np.ndarray:
        """Returns the molecule's Morgan fingerprint"""
        fp_vec = AllChem.GetMorganFingerprintAsBitVect(
            rd_mol,
            self._fingerprint_radius,
            nBits=self._fingerprint_length,
            useChirality=True,
        )
        return np.array(
            fp_vec,
            dtype=bool,
        )
