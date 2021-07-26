""" Module containing classes that implements different cost calculators
"""
from __future__ import annotations
import abc
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem import Molecule


class MoleculeCostCalculator(abc.ABC):
    """
    A base class for all molecule cost calculators.

    The cost can be computed by either calling the `calculate` method
    of by calling the instantiated class with a molecule.

    .. code-block::

        calculator = MyCost()
        cost = calculator.calculate(molecule)
        cost = calculator(molecule)

    :param config: the configuration of the tree search
    """

    def __init__(self, config: Configuration = None) -> None:
        self._config = config

    def __call__(self, mol: Molecule) -> float:
        return self.calculate(mol)

    @abc.abstractmethod
    def calculate(self, mol: Molecule) -> float:
        """
        Calculate the cost of a molecule

        :param mol: the query molecule
        :return: the estimated cost
        """


class ZeroMoleculeCost(MoleculeCostCalculator):
    """Encapsulation of a Zero cost model"""

    def __repr__(self) -> str:
        return "zero"

    def calculate(self, mol: Molecule) -> float:
        return 0.0
