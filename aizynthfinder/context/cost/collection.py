""" Module containing classes to compute the molecule cost in some search algorithms
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.utils.loading import load_dynamic_class
from aizynthfinder.utils.exceptions import CostException
from aizynthfinder.context.cost.costs import MoleculeCostCalculator, ZeroMoleculeCost
from aizynthfinder.context.cost.costs import (
    __name__ as costs_module,
)

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Dict, Any
    from aizynthfinder.chem import Molecule


class MoleculeCost(ContextCollection):
    """ Collection of molecular cost objects """

    _single_selection: bool = True
    _collection_name: str = "molecule cost"

    _aliases: Dict[str, Any] = {
        "zero": ZeroMoleculeCost,
    }

    def __init__(self) -> None:
        super().__init__()
        self._cache: Dict[str, float] = {}
        self._items["zero"] = ZeroMoleculeCost()
        self._selection = ["zero"]

    def __call__(self, mol: Molecule) -> float:
        if not isinstance(self.selection, str):
            raise CostException("No cost selected cannot compute it!")
        if mol.inchi_key not in self._cache:
            self._cache[mol.inchi_key] = self[self.selection](mol)
        return self._cache[mol.inchi_key]

    def load(self, cost: MoleculeCost) -> None:  # type: ignore
        """
        Add a pre-initialized cost calculator object to the collection

        :param cost: the item to add
        """
        if not isinstance(cost, MoleculeCostCalculator):
            raise CostException(
                "Only objects of classes inherited from MoleculeCostCalculator can be added"
            )
        self._items[repr(cost)] = cost

    def load_from_config(self, **costs_config: Any) -> None:
        """
        Load one or more cost models from a configuration

        The format should be:
            key: config_for_model

        :param costs_config: the configuration
        """
        for name_spec, scorer_config in costs_config.items():
            if name_spec in self._aliases:
                cls = self._aliases[name_spec]
            else:
                cls = load_dynamic_class(name_spec, costs_module, CostException)
            obj = cls(**(scorer_config or {}))
            config_str = (
                f" from configuration '{scorer_config}'" if scorer_config else ""
            )
            self._logger.info(f"Loaded cost: '{repr(obj)}'{config_str}")
            self._items[repr(obj)] = obj
