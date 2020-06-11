""" Module containing classes that interfaces different stocks and query classes
"""
import os

import pandas as pd
from pymongo import MongoClient

from aizynthfinder.chem import Molecule
from aizynthfinder.utils.logging import logger


class StockException(Exception):
    """ An exception raised by the Stock classes
    """


class InMemoryInchiKeyQuery:
    """
    A stock query class that is based on an in-memory list
    of pre-computed inchi-keys.

    This list can be instantiated from
        * A Pandas dataframe in HDF5 or CSV format
        * A text file with an inchi key on each row

    The pandas dataframe must have a column "inchi_key" with InChIkeys.
    The HDF5 file must have a dataset called "table".

    :parameter filename: the path to the file with inchi-keys
    :type filename: str
    """

    def __init__(self, filename):
        ext = os.path.splitext(filename)[1]
        if ext in [".h5", ".hdf5"]:
            stock = pd.read_hdf(filename, key="table")
            inchis = stock.inchi_key.values
        elif ext == ".csv":
            stock = pd.read_csv(filename)
            inchis = stock.inchi_key.values
        else:
            with open(filename, "r") as fileobj:
                inchis = [line.strip() for line in fileobj]
        self._stock_inchikeys = set(inchis)

    def __contains__(self, mol):
        return mol.inchi_key in self._stock_inchikeys

    def __len__(self):
        return len(self._stock_inchikeys)

    @property
    def stock_inchikeys(self):
        return self._stock_inchikeys


class MongoDbInchiKeyQuery:
    """
    A stock query class that is looking up inchi keys in a Mongo database.

    The documents in the database collection should have at least 2 fields:
        * inchi_key: the inchi key of the molecule
        * source: the original source of the molecule

    :ivar client: the Mongo client
    :vartype client: pymongo.MongoClient
    :ivar database: the database instance
    :vartype database: pymongo.database.Database
    :ivar molecules: the collection of documents
    :vartype molecules: pymongo.collection.Collection

    :parameter host: the database host, defaults to None
    :type host: str, optional
    :parameter database: the database name, defaults to "stock_db"
    :type database: str, optional
    :parameter collection: the database collection, defaults to "molecules"
    :type collection: str, optional
    """

    def __init__(self, host=None, database="stock_db", collection="molecules"):
        self.client = MongoClient(host or os.environ.get("MONGODB_HOST"))
        self.database = self.client[database]
        self.molecules = self.database[collection]
        self._len = None

    def __contains__(self, mol):
        return self.molecules.count_documents({"inchi_key": mol.inchi_key}) > 0

    def __len__(self):
        if self._len is None:
            self._len = self.molecules.count_documents({})
        return self._len

    def __str__(self):
        return "'MongoDB stock'"

    def availability_string(self, mol):
        """
        Returns the sources of the molecule

        :param mol: the query molecule
        :type mol: Molecule
        :return: a comma-separated list of sources
        :rtype: str
        """
        sources = [
            item["source"] for item in self.molecules.find({"inchi_key": mol.inchi_key})
        ]
        return ",".join(sources)


class Stock:
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

    def __init__(self):
        self._stocks = {}
        self._exclude = set()
        self._selected_stocks = []
        self._logger = logger()

    def __contains__(self, mol):
        if mol.inchi_key in self._exclude:
            return False

        for key in self._selected_stocks:
            if mol in self[key]:
                return True
        return False

    def __getitem__(self, key):
        return self._stocks[key]

    def __len__(self):
        return sum(len(self[key]) for key in self._selected_stocks)

    @property
    def selected_stocks(self):
        """
        Returns the selected sub stocks.
        If set it will update the selection of stocks.

        :return: the selected sub stocks
        :rtype: list of str
        """
        return self._selected_stocks

    @selected_stocks.setter
    def selected_stocks(self, keys):
        self.select_stocks(keys)

    def add_stock(self, key):
        """
        Add a single key to the list of selected stocks

        :param key: the key of the stock
        :type key: str
        """
        if key not in self._stocks:
            raise StockException(f"The {key} stock is not loaded")
        self._selected_stocks.append(key)

    def available_stocks(self):
        """
        Return a list of available stock keys

        :return: the list
        :rtype: list of str
        """
        return list(self._stocks.keys())

    def availability_string(self, mol):
        """
        Return a string of what stocks a given mol is available

        If the molecule is not in stock it will return "Not in stock"

        :param mol: The molecule to query
        :type mol: Molecule

        :returns: string with a list of stocks that mol was found in
        :rtype: str
        """
        availability = []
        for key in self._selected_stocks:
            if mol not in self[key]:
                continue
            if hasattr(self[key], "availability_string"):
                availability.append(self[key].availability_string(mol))
            else:
                availability.append(key)
        if availability:
            return ",".join(availability)
        return "Not in stock"

    def exclude(self, mol):
        """
        Exclude a molecule from the stock.
        When this molecule is queried it will return False,
        regardless if the molecule is in the stock.

        :param mol: the molecule to exclude
        :type mol: Molecule
        """
        self._exclude.add(mol.inchi_key)

    def load_stock(self, source, key):
        """
        Load a stock.

        If `source` is a string, it is taken as a path to a filename and the
        stock is loaded as an `InMemoryInchiKeyQuery` object.

        If `source` is not a string, it is taken as a custom object that
        implements the `__contains__` and `__len__` methods for querying.

        :param source: the source of the sock
        :type source: str or object
        :param key: The key that will be used to select the stock
        :type key: str
        """
        src_str = str(source)
        if "object at 0x" in src_str:
            src_str = source.__class__.__name__
        self._logger.info(f"Loading stock from {src_str} to {key}")

        if isinstance(source, str):
            source = InMemoryInchiKeyQuery(source)
        self._stocks[key] = source

    def remove_stock(self, key):
        """
        Remove a single key from the list of selected stocks

        :param key: The key of the stock
        :type key: str
        """
        if key not in self._selected_stocks:
            raise StockException(
                f"Could not find key {key} in selected stocks: {self.available_stocks}"
            )
        self._selected_stocks.remove(key)

    def reset_exclusion_list(self):
        """Remove all molecules in the exclusion list
        """
        self._exclude = set()

    def select_stocks(self, keys):
        """
        Select what stocks to include by providing a list of keys, or a single key

        :param keys: keys to select
        :type keys: list of str or str
        """
        if isinstance(keys, str):
            keys = [keys]

        for key in keys:
            if key not in self._stocks:
                raise StockException(
                    f"Invalid key specified {key} when selecting stocks"
                )

        self._selected_stocks = list(keys)
        self._logger.info(f"Selected stocks: {', '.join(keys)}")
        try:
            self._logger.info(f"Compounds in stock: {len(self)}")
        except (TypeError, ValueError):  # In case len is not possible to compute
            pass

    def smiles_in_stock(self, smiles):
        """
        Check if the SMILES is in the currently selected stocks

        :param smiles: SMILES string (must be RDKit sanitizable)
        :type smiles: str
        :returns: if the SMILES was on stock
        :rtype: bool
        """
        return Molecule(smiles=smiles) in self
