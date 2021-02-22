Command-line interface
======================

This tools provide the possibility to perform tree search on a batch of molecules.

In its simplest form, you type

.. code-block:: bash

    aizynthcli --config config_local.yml --smiles smiles.txt

where `config_local.yml` contains configurations such as paths to policy models and stocks (see :doc:`here <configuration>`)
and `smiles.txt` is a simple text file with SMILES (one on each row). 


To find out what other arguments are available use the ``-h`` flag.

.. code-block:: bash

    aizynthcli -h

That gives something like this:

.. include:: cli_help.txt


Analysing output
----------------


The results from the ``aizynthcli`` tool when supplying multiple SMILES is an HDF5 file that can be read as a pandas dataframe. 
It will be called `output.hfdf5` by default.

.. code-block::

  import pandas as pd
  data = pd.read_hdf("output.hdf5", "table")

it will contain statistics about the tree search and the top-ranked routes (as JSONs) for each target compund.

When a single SMILES is provided to the tool, the statistics will be written to the terminal, and the top-ranked routes to
a JSON file (`trees.json` by default).


This is an example of how to create images of the top-ranked routes for the first target compound


.. code-block::

    import pandas as pd
    from aizynthfinder.analysis import ReactionTree

    data = pd.read_hdf("output.hdf5", "table")
    all_trees = data.trees.values  # This contains a list of all the trees for all the compounds
    trees_for_first_target = all_trees[0]

    for itree, tree in enumerate(trees_for_first_target):
        imagefile = f"route{itree:03d}.png"
        ReactionTree.from_dict(tree).to_image().save(imagefile)

The images will be called `route000.png`, `route001.png` etc. 