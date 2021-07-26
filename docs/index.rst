aizynthfinder documentation
===========================

aizynthfinder is a tool for retrosynthetic planning. The algorithm is based on a Monte Carlo tree search that recursively breaks down a molecule to purchasable precursors. The tree search is guided by a policy that suggests possible precursors by utilizing a neural network trained on a library of known reaction templates.  

Introduction
------------

There are two main interfaces provided by the package:
a script that performs tree search in batch mode and
an interface that is providing a GUI within a Jupyter notebook.


The GUI interface should be run in a Jupyter notebook. This is a simple example of the code in a Jupyter notebook cell.

.. code-block::

    from aizynthfinder.interfaces import AiZynthApp
    app = AiZynthApp("/path/to/configfile.yaml")

where the ``AiZynthApp`` class needs to be instantiated with the path to a configuration file (see :doc:`here <configuration>`).

To use the interface, follow these steps:

  1. Executed the code in the cell (press ``Ctrl+Enter``) and a simple GUI will appear
  2. Enter the target SMILES and select stocks and policy model. 
  3. Press the ``Run Search`` button to perform the tree search.
  4. Press the ``Show Reactions`` to see the top-ranked routes



The batch-mode script is called ``aizynthcli`` and can be executed like:

.. code-block:: bash

    aizynthcli --config config_local.yml --policy full_uspto --stock ACD eMolecules --smiles smiles.txt


where `config_local.yml` contains configurations such as paths to policy models and stocks (see :doc:`here <configuration>`),
`full_uspto` is the chosen policy model, `ACD` and `eMolecules` are the chosen stocks and `smiles.txt` is a simple text
file with SMILES (one on each row). If the ``--policy`` or ``--stock`` arguments are omitted, the first policy and all stocks
in the configuration file are auto-selected.

If you just want to perform the tree search on a single molecule. You can directly specify it on the command-line
within quotes:

.. code-block:: bash

    aizynthcli --config config_local.yml --policy full_uspto --stock ACD Addendum --smiles "COc1cccc(OC(=O)/C=C/c2cc(OC)c(OC)c(OC)c2)c1"


The output is some statistics about the tree search, the scores of the top-ranked routes, and the reaction tree
of the top-ranked routes. When smiles are provided in a text file the results are stored in an HDF5 file,
whereas if the SMILEs is provided on the command-line it is printed directly to the prompt 
(except the reaction trees, which are written to a JSON file).


.. toctree::
    :hidden:
    
    gui
    cli
    python_interface
    configuration
    stocks
    scoring
    training    
    aizynthfinder
    sequences
    relationships