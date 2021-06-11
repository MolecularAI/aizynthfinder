""" Module containing classes to compute the molecule cost in some search algorithms
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.context.collection import ContextCollection

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Dict, Any
    from aizynthfinder.chem import Molecule


class CostException(Exception):
    """Exception raised by classes in this module"""


class ZeroCost:
    """Encapsulation of a Zero cost model"""

    def __call__(self, _) -> float:
        return 0.0

    def __repr__(self) -> str:
        return "zero"


class MoleculeCost(ContextCollection):
    """ Collection of molecular cost objects """
    _single_selection: bool = True
    _collection_name: str = "molecule cost"

    _recognized_costs: Dict[str, Any] = {
        "zero": ZeroCost,
    }

    def __init__(self) -> None:
        super().__init__()
        self._cache: Dict[str, float] = {}
        self._items["zero"] = ZeroCost()
        self._selection = ["zero"]

    def __call__(self, mol: Molecule) -> float:
        if not isinstance(self.selection, str):
            raise CostException("No cost selected cannot compute it!")
        if mol.inchi_key not in self._cache:
            self._cache[mol.inchi_key] = self[self.selection](mol)
        return self._cache[mol.inchi_key]

    def load(self, cost: Any) -> None:  # type: ignore
        """
        Load a molecule cost under the given key

        :param cost: A molecule cost object that can be called
        """
        try:
            _ = cost(None)
        except AttributeError:
            raise CostException(
                "Only objects of classes that return a float when caled can be used as costs"
            )
        except Exception:  # pylint: disable=broad-except
            pass
        self._items[repr(cost)] = cost

    def load_from_config(self, **costs_config: Any) -> None:
        """
        Load one or more cost models from a configuration

        The format should be:
            key: config_for_model

        :param config: the configuration
        :type config: key value pairs
        """
        for name_spec, scorer_config in costs_config.items():
            config_str = (
                f" from configuration '{scorer_config}'" if scorer_config else ""
            )
            if name_spec in self._recognized_costs:
                cls = self._recognized_costs[name_spec]
            else:
                cls = self._load_dynamic_cls(name_spec, self.__module__, CostException)
            obj = cls(**(scorer_config or {}))
            self._logger.info(f"Loaded cost: '{repr(obj)}'{config_str}")
            self._items[repr(obj)] = obj
