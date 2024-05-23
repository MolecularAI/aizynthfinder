""" Module containing Retro* cost models """
from __future__ import annotations

import pickle
from typing import TYPE_CHECKING

import numpy as np

from aizynthfinder.search.retrostar.cost import __name__ as retrostar_cost_module
from aizynthfinder.utils.loading import load_dynamic_class

if TYPE_CHECKING:
    from aizynthfinder.chem import Molecule
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Any, List, Tuple


class MoleculeCost:
    """
    A class to compute the molecule cost.

    The cost to be computed is taken from the input config. If no `molecule_cost` is
    set, assigns ZeroMoleculeCost as the `cost` by default. The `molecule_cost` can be
    set as a dictionary in config under `search` in the following format:
    'algorithm': 'retrostar'
    'algorithm_config': {
        'molecule_cost': {
            'cost': name of the search cost class or custom_package.custom_model.CustomClass,
            other settings or params
        }
    }

    The cost can be computed by calling the instantiated class with a molecule.

    .. code-block::

        calculator = MyCost(config)
        cost = calculator.calculate(molecule)

    :param config: the configuration of the tree search
    """

    def __init__(self, config: Configuration) -> None:
        self._config = config
        if "molecule_cost" not in self._config.search.algorithm_config:
            self._config.search.algorithm_config["molecule_cost"] = {
                "cost": "ZeroMoleculeCost"
            }
        kwargs = self._config.search.algorithm_config["molecule_cost"].copy()

        cls = load_dynamic_class(kwargs["cost"], retrostar_cost_module)
        del kwargs["cost"]

        self.molecule_cost = cls(**kwargs) if kwargs else cls()

    def __call__(self, mol: Molecule) -> float:
        return self.molecule_cost.calculate(mol)


class RetroStarCost:
    """
    Encapsulation of the original Retro* molecular cost model

    Numpy implementation of original pytorch model

    The predictions of the score is made on a Molecule object

    .. code-block::

        mol = Molecule(smiles="CCC")
        scorer = RetroStarCost()
        score = scorer.calculate(mol)

    The model provided when creating the scorer object should be a pickled
    tuple.
    The first item of the tuple should be a list of the model weights for each layer.
    The second item of the tuple should be a list of the model biases for each layer.

    :param model_path: the filename of the model weights and biases
    :param fingerprint_length: the number of bits in the fingerprint
    :param fingerprint_radius: the radius of the fingerprint
    :param dropout_rate: the dropout_rate
    """

    _required_kwargs = ["model_path"]

    def __init__(self, **kwargs: Any) -> None:
        model_path = kwargs["model_path"]
        self.fingerprint_length: int = int(kwargs.get("fingerprint_length", 2048))
        self.fingerprint_radius: int = int(kwargs.get("fingerprint_radius", 2))
        self.dropout_rate: float = float(kwargs.get("dropout_rate", 0.1))

        self._dropout_prob = 1.0 - self.dropout_rate
        self._weights, self._biases = self._load_model(model_path)

    def __repr__(self) -> str:
        return "retrostar"

    def calculate(self, mol: Molecule) -> float:
        # pylint: disable=invalid-name
        mol.sanitize()
        vec = mol.fingerprint(
            radius=self.fingerprint_radius, nbits=self.fingerprint_length
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

    @staticmethod
    def _load_model(model_path: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        with open(model_path, "rb") as fileobj:
            weights, biases = pickle.load(fileobj)

        return (
            [np.asarray(item) for item in weights],
            [np.asarray(item) for item in biases],
        )


class ZeroMoleculeCost:
    """Encapsulation of a Zero cost model"""

    def __repr__(self) -> str:
        return "zero"

    def calculate(self, _mol: Molecule) -> float:  # pytest: disable=unused-argument
        return 0.0
