""" Module containing classes that interfaces different stocks query classes
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import pandas as pd

try:
    import molbloom
except ImportError:
    HAS_MOLBLOOM = False
else:
    HAS_MOLBLOOM = True

from aizynthfinder.chem import Molecule
from aizynthfinder.utils.exceptions import StockException
from aizynthfinder.utils.mongo import get_mongo_client

if TYPE_CHECKING:
    from pymongo.collection import Collection as MongoCollection
    from pymongo.database import Database as MongoDatabase

    from aizynthfinder.utils.type_utils import Optional, Set, StrDict


class StockQueryMixin:
    """
    Mixin class for all query classes, providing a default interface
    to some methods that might not be possible to implement for each
    query class.
    """

    def __len__(self) -> int:
        return 0

    def __contains__(self, mol: Molecule) -> bool:
        return False

    def amount(self, mol: Molecule) -> float:
        """
        Returns the maximum amount of the molecule in stock

        :param mol: the query molecule
        :raises StockException: if the amount cannot be computed
        :return: the amount
        """
        raise StockException("Cannot compute amount")

    def availability_string(self, mol: Molecule) -> str:
        """
        Returns the sources of the molecule

        :param mol: the query molecule
        :raises StockException: if the string cannot be computed
        :return: a comma-separated list of sources
        """
        raise StockException("Cannot provide availability")

    def cached_search(self, mol: Molecule) -> bool:
        """
        Finds the entries of the molecule in the stock and cache them
        if necessary.

        :param mol: the query molecule
        :return: if the molecule is in stock
        """
        return mol in self

    def clear_cache(self) -> None:
        """Clear the internal search cache if available"""

    def price(self, mol: Molecule) -> float:
        """
        Returns the minimum price of the molecule in stock

        :param mol: the query molecule
        :raises StockException: if the price cannot be computed
        :rtype: float
        """
        raise StockException("Cannot compute price")


class InMemoryInchiKeyQuery(StockQueryMixin):
    """
    A stock query class that is based on an in-memory list
    of pre-computed inchi-keys.

    This list can be instantiated from
        * A Pandas dataframe in HDF5 or CSV format
        * A text file with an inchi key on each row

    The dataframe must have a column with InChIkeys that by default is "inchi_key".
    The HDF5 file must have a dataset called "table".

    If the source is a dataframe, then optionally it can contain prices and this
    columns can be specified with the "price_column" argument.

    :parameter path: the path to the file with inchi-keys
    :parameter inchi_key_col: the name of the column of the InChI keys
    :paramater price_col: the name of the column with the optional prices
    """

    def __init__(
        self,
        path: str,
        inchi_key_col: str = "inchi_key",
        price_col: Optional[str] = None,
    ) -> None:
        ext = os.path.splitext(path)[1]
        if ext not in [".h5", ".hdf5", ".csv"]:
            with open(path, "r") as fileobj:
                inchis = fileobj.read().splitlines()
            self._stock_inchikeys = set(inchis)
            self._price_dict: StrDict = {}
            return

        if ext in [".h5", ".hdf5"]:
            stock_df: pd.DataFrame = pd.read_hdf(path, key="table")
        else:
            stock_df = pd.read_csv(
                path,
                usecols=[inchi_key_col, price_col] if price_col else [inchi_key_col],
            )
        inchis = stock_df[inchi_key_col].values
        self._stock_inchikeys = set(inchis)

        if price_col is None:
            self._price_dict = {}
            return

        if len(stock_df) != len(self._stock_inchikeys):
            raise StockException(
                "InChI keys in stock df are expected to be unique, currently they are not"
            )
        if stock_df[price_col].isnull().sum() != 0:
            raise StockException(
                "null values found in the price column, please drop/impute them first"
            )
        if stock_df[price_col].min() < 0:
            raise StockException(
                f"expected non-negative prices, the min in the current file {path} is {stock_df[price_col].min()} "
            )
        self._price_dict = dict(zip(stock_df[inchi_key_col], stock_df[price_col]))

    def __contains__(self, mol: Molecule) -> bool:
        return mol.inchi_key in self._stock_inchikeys

    def __len__(self) -> int:
        return len(self._stock_inchikeys)

    @property
    def stock_inchikeys(self) -> Set[str]:
        """Return the InChiKeys in this stock"""
        return self._stock_inchikeys

    def price(self, mol: Molecule) -> float:
        if not self._price_dict:
            raise StockException(
                "no prices created, check the path type and if price column is supplied"
            )
        if mol in self:
            return self._price_dict[mol.inchi_key]
        raise StockException(f"no price info available for {mol.smiles}")


class MongoDbInchiKeyQuery(StockQueryMixin):
    """
    A stock query class that is looking up inchi keys in a Mongo database.

    The documents in the database collection should have at least 2 fields:
        * inchi_key: the inchi key of the molecule
        * source: the original source of the molecule

    :ivar client: the Mongo client
    :ivar database: the database instance
    :ivar molecules: the collection of documents

    :parameter host: the database host, defaults to None
    :parameter database: the database name, defaults to "stock_db"
    :parameter collection: the database collection, defaults to "molecules"
    """

    def __init__(
        self,
        host: Optional[str] = None,
        database: str = "stock_db",
        collection: str = "molecules",
    ) -> None:
        self.client = get_mongo_client(
            host or os.environ.get("MONGODB_HOST") or "localhost"
        )
        if self.client is None:
            raise ImportError(
                "Cannot use this stock query class because it seems like pymongo is not installed. "
                "Please install aizynthfinder with extras dependencies."
            )
        self.database: MongoDatabase = self.client[database]
        self.molecules: MongoCollection = self.database[collection]
        self._len: Optional[int] = None

    def __contains__(self, mol: Molecule) -> bool:
        return self.molecules.count_documents({"inchi_key": mol.inchi_key}) > 0

    def __len__(self) -> int:
        if self._len is None:
            self._len = self.molecules.count_documents({})
        return self._len

    def __str__(self) -> str:
        return "'MongoDB stock'"

    def availability_string(self, mol: Molecule) -> str:
        sources = [
            item["source"] for item in self.molecules.find({"inchi_key": mol.inchi_key})
        ]
        return ",".join(sources)


class MolbloomFilterQuery(StockQueryMixin):
    """
    A stock query class that is based on an a molbloom filter
    for SMILES strings or InChI keys

    :parameter path: the path to the saved bloom filter
    :parameter smiles_based: if True will use SMILES for lookup instead of InChI keys
    """

    def __init__(self, path: str, smiles_based: bool = False) -> None:
        if not HAS_MOLBLOOM:
            raise ImportError(
                "Cannot use this stock query class because it seems like molbloom is not installed. "
                "Please install aizynthfinder with extras dependencies."
            )
        self._filter = molbloom.BloomFilter(path)
        self._smiles_based = smiles_based

    def __contains__(self, mol: Molecule) -> bool:
        if self._smiles_based:
            return mol.smiles in self._filter
        return mol.inchi_key in self._filter


STOCK_QUERY_ALIAS = {
    "inchiset": "InMemoryInchiKeyQuery",
    "mongodb": "MongoDbInchiKeyQuery",
    "bloom": "MolbloomFilterQuery",
}
