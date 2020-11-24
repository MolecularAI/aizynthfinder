Python interface
================

This page gives a quick example of how the tree search can be completed
by writing your own python interface. This is not recommended for most users.


1. Import the necessary class

.. code-block:: python

    from aizynthfinder.aizynthfinder import AiZynthFinder


2. Instantiate that class by providing a configuration file.


.. code-block:: python

    filename = "config.yml"
    finder = AiZynthFinder(configfile=filename)


3. Select stock and policy


.. code-block:: python

    finder.stock.select("zinc")
    finder.policy.select("full_uspto")

`zinc` and `full_uspto` where the keys given to the stock and the policy in the configuration file.


4. Set the target SMILES and perform the tree search


.. code-block:: python

    finder.target_smiles = "Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C"
    finder.tree_search()


5. Analyse the search tree and build routes


.. code-block:: python

    finder.build_routes()
    stats = finder.extract_statistics()


The ``build_routes`` method needs to be called before any analysis can be done.


Further reading
---------------

The docstrings of all modules, classes and methods can be consulted :doc:`here <aizynthfinder>`


and you can always find them in an interactive Python shell using for instance:

.. code-block:: python

    from aizynthfinder.chem import Molecule
    help(Molecule)
    help(Molecule.fingerprint)


If you are interested in the the relationships between the classes have a look :doc:`here <relationships>`
and if you want to dig deeper in to the main algorithmic sequences have a look :doc:`here <sequences>`