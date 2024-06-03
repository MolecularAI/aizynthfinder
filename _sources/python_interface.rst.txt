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
    finder.expansion_policy.select("uspto")
    finder.filter_policy.select("uspto")

`zinc` and `uspto` where the keys given to the stock and the policy in the configuration file.
The first policy set is the expansion policy and the second is the filter policy. The filter policy is optional.

4. Set the target SMILES and perform the tree search


.. code-block:: python

    finder.target_smiles = "Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C"
    finder.tree_search()


5. Analyse the search tree and build routes


.. code-block:: python

    finder.build_routes()
    stats = finder.extract_statistics()


The ``build_routes`` method needs to be called before any analysis can be done.

Expansion interface
-------------------

There is an interface for the expansion policy as well. It can be used to break down a molecule into reactants.

.. code-block:: python

    filename = "config.yml"
    expander = AiZynthExpander(configfile=filename)
    expander.expansion_policy.select("uspto")
    expander.filter_policy.select("uspto")
    reactions = expander.do_expansion("Cc1cccc(c1N(CC(=O)Nc2ccc(cc2)c3ncon3)C(=O)C4CCS(=O)(=O)CC4)C")

for this, you only need to select the policies. The filter policy is optional and using it will only add the
feasibility of the reactions not filter it out.

The result is a nested list of `FixedRetroReaction` objects. This you can manipulate to for instance get
out all the reactants SMILES strings

.. code-block:: python

    reactants_smiles = []
    for reaction_tuple in reactions:
        reactants_smiles.append([mol.smiles for mol in reaction_tuple[0].reactants[0])

or you can put all the metadata of all the reactions in a pandas dataframe

.. code-block:: python

    import pandas as pd
    metadata = []
    for reaction_tuple in reactions:
        for reaction in reaction_tuple:
            metadata.append(reaction.metadata)
    df = pd.DataFrame(metadata)


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