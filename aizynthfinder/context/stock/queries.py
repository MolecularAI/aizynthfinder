""" Module containing classes that interfaces different stocks query classes
"""
from __future__ import annotations
import os
from typing import TYPE_CHECKING

import pandas as pd

from aizynthfinder.chem import Molecule
from aizynthfinder.utils.mongo import get_mongo_client
from aizynthfinder.utils.exceptions import StockException

if TYPE_CHECKING:
    from pymongo.database import Database as MongoDatabase
    from pymongo.collection import Collection as MongoCollection

    # pylint: disable=ungrouped-imports
    from aizynthfinder.utils.type_utils import Set, Optional


# pylint: disable=no-self-use
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

    The pandas dataframe must have a column "inchi_key" with InChIkeys.
    The HDF5 file must have a dataset called "table".

    :parameter filename: the path to the file with inchi-keys
    """

    def __init__(self, filename: str) -> None:
        ext = os.path.splitext(filename)[1]
        if ext in [".h5", ".hdf5"]:
            stock = pd.read_hdf(filename, key="table")  # type: ignore
            inchis = stock.inchi_key.values  # type: ignore
        elif ext == ".csv":
            stock = pd.read_csv(filename)
            inchis = stock.inchi_key.values
        else:
            with open(filename, "r") as fileobj:
                inchis = fileobj.read().splitlines()
        self._stock_inchikeys = set(inchis)

    def __contains__(self, mol: Molecule) -> bool:
        return mol.inchi_key in self._stock_inchikeys

    def __len__(self) -> int:
        return len(self._stock_inchikeys)

    @property
    def stock_inchikeys(self) -> Set[str]:
        """Return the InChiKeys in this stock"""
        return self._stock_inchikeys


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
        host: str = None,
        database: str = "stock_db",
        collection: str = "molecules",
    ) -> None:
        self.client = get_mongo_client(
            host or os.environ.get("MONGODB_HOST") or "localhost"
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
