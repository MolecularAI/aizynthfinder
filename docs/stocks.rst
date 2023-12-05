Stocks
======

The stock files specified in the configuration file are loaded and a set of inchi keys
are stored in-memory for lookup. However, the tool supports other stock queries as well as a way
to fully customize the lookup.

Mongo database stock
--------------------

First, support for lookup inchi keys in a Mongo database is supported. The Mongo client should
have a database and a collection containing documents with at least two fields: `inchi_key` and `source`.
The `inchi_key` field will be used for lookup and `source` specifies the source database of the compound.

By adding these lines to the configuration file, the Mongo database will be used:

.. code-block:: yaml

    stock:
        type: mongodb
        host: user@myurl.com
        database: database_name
        collection: compounds


If no options are provided to the ``mongodb_stock`` key, the host, database and collection are taken to be `localhost`, 
`stock_db`, and `molecules`, respectively. 

Stop criteria
-------------

The stock can be used to stop the tree search based on three criteria: a) minimum price, b) maximum amount and c) count of different elements in the molecule.
Note that the stock query class need to support querying for price and amount, if the stop criteria should work properly.

The stop criteria can be specified in the configuration file 

.. code-block:: yaml

    stock:
        stop_criteria:
            price: 10
            counts:
                C: 10


In the Jupyter GUI you can set the limit on the element occurences, but currently not the price and amount limits. 

Custom stock
------------

Support for any type of lookup is provided. You just need to write a python class that implements the ``__contains__`` 
and subclasses the ``aizynthfinder.context.stock.queries.StockQueryMixin``. The ``__contains__`` method is used for lookup and should take a ``Molecule`` object as only argument.
The ``StockQueryMixin`` mixin class provide a default interface for some methods that perhaps isn't possible to implement in all query classes.

This is an example:

.. code-block::

  from rdkit.Chem import Lipinski
  from aizynthfinder.context.stock.queries import StockQueryMixin
  class CriteriaStock(StockQueryMixin):
      def __contains__(self, mol):
          return Lipinski.HeavyAtomCount(mol.rd_mol) < 10


To use this stock with the ``aizynthcli`` tool, save it in a ``custom_stock.py`` module that is located in a directory known to 
the python interpreter. Add this line to the module.

.. code-block::

  stock = CriteriaStock()


and it will be automatically used in the tree search. 

Alternatively the custom query class can be used by the ``aizynthapp`` tool.


.. code-block::

  from aizynthfinder import AiZynthApp
  configfile="config_local.yml"
  app = AiZynthApp(configfile, setup=False)
  app.finder.stock.load(CriteriaStock(), "criteria") # This loads the custom stock class
  app.setup()


Lastly, it is possible to specify a custom stock class in the configuration file if it is located in a module that 
is known by the python interpreter.

.. code-block::

    stock:
        type: aizynthfinder.contrib.stocks.CriteriaStock


can be used if the `aizynthfinder.contrib.stocks` is an existing sub-package and module.


Making stocks
-------------

We provide a tool to create inchi key-based stocks from SMILES strings. Thereby, one
can create a stock based on for instance a subset of the ZINC database.

The tool support both creating a stock in HDF5 format or adding them to an existing Mongo database.

The tool is easiest to use if one has a number of plain text files, in which each row has one SMILES.

Then one can use one of these two commands:


.. code-block::

    smiles2stock --files file1.smi file2.smi --output stock.hdf5
    smiles2stock --files file1.smi file2.smi --output my_db --target mongo


to create either an HDF5 stock or a Mongo database stock, respectively. The ``file1.smi`` and ``file2.smi``
are simple text files and ``my_db`` is the source tag for the Mongo database.


If one has SMILES in any other format, one has to provide a custom module that extract the SMILES from
the input files. This is an example of such a module that can be used with downloads from the Zinc database
where the first row contains headers and the SMILES are the first element on each line.


.. code-block::

    def extract_smiles(filename):
        with open(filename, "r") as fileobj:
            for i, line in enumerate(fileobj.readlines()):
                if i == 0:
                    continue
                yield line.strip().split(" ")[0]


if this is saved as ``load_zinc.py`` in a path that is known to the Python interpreter, it can be 
used like this

.. code-block::

    export PYTHONPATH=`pwd`
    smiles2stock --files load_zinc file1.smi file2.smi --source module --output stock.hdf5


where the first line adds the current directory to the python path (if you are using a Bash shell).
