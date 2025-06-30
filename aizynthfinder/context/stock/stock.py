""" Module containing classes that interfaces different stock classes
"""

from __future__ import annotations

import copy
from collections import defaultdict
from typing import TYPE_CHECKING

from aizynthfinder.chem import Molecule
from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.context.stock.queries import (
    STOCK_QUERY_ALIAS,
    InMemoryInchiKeyQuery,
    MolbloomFilterQuery,
    StockQueryMixin,
)
from aizynthfinder.context.stock.queries import __name__ as queries_module
from aizynthfinder.utils.exceptions import StockException
from aizynthfinder.utils.loading import load_dynamic_class

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        List,
        Optional,
        Set,
        StrDict,
        Union,
    )


class Stock(ContextCollection):
    """
    A collection of molecules that are in stock

    A molecule can be queried on the stock with:

    .. code-block::

        my_mol = Molecule(smiles="CCO")
        my_mol in stock

    One can obtain individual stocks with:

    .. code-block::

        sub_stock = stock["key"]

    One can obtain the number of molecules in the selected stock with:

    .. code-block::

        number_of_molecules = len(stock)

    """

    _collection_name = "stock"

    def __init__(self) -> None:
        super().__init__()
        self._exclude: Set[str] = set()
        self._stop_criteria: StrDict = {
            "amount": None,
            "price": None,
            "weight": None,
            "counts": {},
        }
        self._use_stop_criteria: bool = False

    def __contains__(self, mol: Molecule) -> bool:
        if not self.selection or mol.inchi_key in self._exclude:
            return False

        if self._use_stop_criteria:
            return self._apply_stop_criteria(mol)

        for key in self.selection:
            if mol in self[key]:
                return True
        return False

    def __len__(self) -> int:
        return sum(len(self[key]) for key in self.selection or [])

    @property
    def stop_criteria(self) -> dict:
        """Return a copy of the stop criteria used by the stock"""
        return copy.deepcopy(self._stop_criteria)

    def amount(self, mol: Molecule) -> float:
        """
        Calculate the maximum amount of a molecule in stock

        :param mol: the molecule to query
        :raises StockException: if the amount could not be computed
        :return: the maximum amount
        """
        amounts = self._mol_property(mol, "amount")
        if not amounts:
            raise StockException("Could not obtain amount of molecule")
        return max(amounts)

    def availability_list(self, mol: Molecule) -> List[str]:
        """
        Return a list of what stocks a given mol is available

        If the molecule is not in stock it will return any empty list

        :param mol: The molecule to query
        :returns: string with a list of stocks that mol was found in
        """
        availability = []
        for key in self.selection or []:
            if mol not in self[key]:
                continue
            try:
                availability.extend(self[key].availability_string(mol).split(","))
            except (StockException, AttributeError):
                availability.append(key)
        return sorted(list(set(availability)))

    def availability_string(self, mol: Molecule) -> str:
        """
        Return a string of what stocks a given mol is available

        If the molecule is not in stock it will return "Not in stock"

        :param mol: The molecule to query
        :returns: string with a list of stocks that mol was found in
        """
        availability = self.availability_list(mol)
        if availability:
            return ",".join(availability)
        return "Not in stock"

    def exclude(self, mol: Molecule) -> None:
        """
        Exclude a molecule from the stock.
        When this molecule is queried it will return False,
        regardless if the molecule is in the stock.

        :param mol: the molecule to exclude
        """
        self._exclude.add(mol.inchi_key)

    def load(self, source: StockQueryMixin, key: str) -> None:  # type: ignore
        """
        Add a pre-initialized stock query object to the stock

        :param source: the item to add
        :param key: The key that will be used to select the stock
        """
        if not isinstance(source, StockQueryMixin):
            raise StockException(
                "Only objects of classes inherited from StockQueryMixin can be added"
            )

        self._logger.info(f"Loading stock from {source.__class__.__name__} to {key}")
        self._items[key] = source

    def load_from_config(self, **config: Any) -> None:
        """
        Load one or more stock queries from a configuration

        The key can be "stop_criteria" in case the config is given to the
        `set_stop_criteria` method

        The format should be
        key:
            type: name of the stock class or custom_package.custom_model.CustomClass
            path: path to the stock file
            other settings or params
        or
        key: path_to_model

        :param config: the configuration
        """
        if "stop_criteria" in config:
            self.set_stop_criteria(config["stop_criteria"])

        for key, stock_config in config.items():
            if key == "stop_criteria":
                continue

            if not isinstance(stock_config, dict):
                kwargs = {"path": stock_config}
                if stock_config.endswith(".bloom"):
                    cls: Any = MolbloomFilterQuery
                else:
                    cls = InMemoryInchiKeyQuery
            else:
                if "type" not in stock_config or stock_config["type"] == "inchiset":
                    cls = InMemoryInchiKeyQuery
                else:
                    stock_query = STOCK_QUERY_ALIAS.get(
                        stock_config["type"], stock_config["type"]
                    )
                    cls = load_dynamic_class(
                        stock_query, queries_module, StockException
                    )
                kwargs = dict(stock_config)

            if "type" in kwargs:
                del kwargs["type"]
            obj = cls(**kwargs)
            self.load(obj, key)

    def price(self, mol: Molecule) -> float:
        """
        Calculate the minimum price of a molecule in stock

        :param mol: the molecule to query
        :raises StockException: if the price could not be computed
        :return: the minimum price
        """
        prices = self._mol_property(mol, "price")
        if not prices:
            raise StockException("Could not obtain price of molecule")
        return min(prices)

    def reset_exclusion_list(self) -> None:
        """Remove all molecules in the exclusion list"""
        self._exclude = set()

    def select(self, value: Union[str, List[str]], append: bool = False) -> None:
        """
        Select one or more stock queries

        :param value: the key of the stocks to select
        :param append: if True and ``value`` is a single key append it to the current selection
        """
        super().select(value, append)
        try:
            self._logger.info(f"Compounds in stock: {len(self)}")
        except (TypeError, ValueError):  # In case len is not possible to compute
            pass

    def set_stop_criteria(self, criteria: Optional[Dict] = None) -> None:
        """
        Set criteria that stop the search

        The keys of the criteria can be "price", "amount", or "weight" which accepts numerical settings,
        or "counts" that should be dictionary of maximum allowed count for each atomic symbol.

        Example:

        .. code-block::

            criteria = {
                "price": 5,
                "amount": 100,
                "weight": 250,
                "counts": {
                    "C": 6,
                    "O": 4
                }
            }

        :param criteria: the criteria settings
        """
        criteria = criteria or {}
        self._stop_criteria = {
            "price": criteria.get("price"),
            "amount": criteria.get("amount"),
            "weight": criteria.get("weight"),
            "counts": copy.deepcopy(criteria.get("size", criteria.get("counts"))),
        }
        self._use_stop_criteria = any(self._stop_criteria.values())
        reduced_criteria = {
            key: value for key, value in self._stop_criteria.items() if value
        }
        self._logger.info(f"Stop criteria for stock updated to: {reduced_criteria}")

    def smiles_in_stock(self, smiles: str) -> bool:
        """
        Check if the SMILES is in the currently selected stocks

        :param smiles: SMILES string (must be RDKit sanitizable)
        :returns: if the SMILES was on stock
        """
        return Molecule(smiles=smiles) in self

    def _apply_amount_criteria(self, mol: Molecule) -> bool:
        if not self._stop_criteria["amount"]:
            return True
        try:
            amount = self.amount(mol)
        except StockException:
            return True
        return amount >= self._stop_criteria.get("amount", amount)

    def _apply_counts_criteria(self, mol: Molecule) -> bool:
        if not self._stop_criteria["counts"]:
            return True
        atom_counts: dict = defaultdict(int)
        for atom in mol.rd_mol.GetAtoms():
            atom_counts[atom.GetSymbol()] += 1
        for symbol, threshold in self._stop_criteria["counts"].items():
            if atom_counts[symbol] > threshold:
                return False
        return True

    def _apply_price_criteria(self, mol: Molecule) -> bool:
        if not self._stop_criteria["price"]:
            return True
        try:
            price = self.price(mol)
        except StockException:
            return True
        return price <= self._stop_criteria.get("price", price)

    def _apply_stop_criteria(self, mol: Molecule) -> bool:
        if not self._apply_weight_criteria(mol):
            return False

        if not self._apply_counts_criteria(mol):
            return False

        passes = False
        for key in self.selection or []:
            passes = passes or self[key].cached_search(mol)
        passes = passes and self._apply_amount_criteria(mol)
        passes = passes and self._apply_price_criteria(mol)

        for key in self.selection or []:
            self[key].clear_cache()
        return passes

    def _apply_weight_criteria(self, mol: Molecule) -> bool:
        if not self._stop_criteria["weight"]:
            return True
        return mol.weight <= self._stop_criteria["weight"]

    def _mol_property(self, mol, property_name):
        values = []
        for key in self.selection:
            try:
                func = getattr(self[key], property_name)
                values.append(func(mol))
            except StockException:
                pass
        return values
