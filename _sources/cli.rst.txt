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


By default:

  * `All` stocks are selected if no stock is specified
  * `First` expansion policy is selected if not expansion policy is specified
  * `All` filter policies are selected if it is not specified on the command-line

Analysing output
----------------


The results from the ``aizynthcli`` tool when supplying multiple SMILES is an JSON or HDF5 file that can be read as a pandas dataframe.
It will be called `output.json.gz` by default.

A `checkpoint.json.gz` will also be generated if a checkpoint file path is provided as input when calling the ``aizynthcli`` tool. The
checkpoint data will contain the processed smiles with their corresponding results in each line of the file.

.. code-block::

  import pandas as pd
  data = pd.read_json("output.json.gz", orient="table")

it will contain statistics about the tree search and the top-ranked routes (as JSONs) for each target compound, see below.

When a single SMILES is provided to the tool, the statistics will be written to the terminal, and the top-ranked routes to
a JSON file (`trees.json` by default).


This is an example of how to create images of the top-ranked routes for the first target compound


.. code-block::

    import pandas as pd
    from aizynthfinder.reactiontree import ReactionTree

    data = pd.read_json("output.json.gz", orient="table")
    all_trees = data.trees.values  # This contains a list of all the trees for all the compounds
    trees_for_first_target = all_trees[0]

    for itree, tree in enumerate(trees_for_first_target):
        imagefile = f"route{itree:03d}.png"
        ReactionTree.from_dict(tree).to_image().save(imagefile)

The images will be called `route000.png`, `route001.png` etc.


Specification of output
-----------------------

The JSON or HDF5 file created when running the tool with a list of SMILES will have the following columns

============================= ===========
Column                        Description
============================= ===========
target                        The target SMILES
search_time                   The total search time in seconds
first_solution_time           The time elapsed until the first solution was found
first_solution_iteration      The number of iterations completed until the first solution was found
number_of_nodes               The number of nodes in the search tree
max_transforms                The maximum number of transformations for all routes in the search tree
max_children                  The maximum number of children for a search node
number_of_routes              The number of routes in the search tree
number_of_solved_routes       The number of solved routes in search tree
top_score                     The score of the top-scored route (default to MCTS reward)
is_solved                     If the top-scored route is solved
number_of_steps               The number of reactions in the top-scored route
number_of_precursors          The number of starting materials
number_of_precursors_in_stock The number of starting materials in stock
precursors_in_stock           Comma-separated list of SMILES of starting material in stock
precursors_not_in_stock       Comma-separated list of SMILES of starting material not in stock
precursors_availability       Semi-colon separated list of stock availability of the staring material
policy_used_counts            Dictionary of the total number of times an expansion policy have been used
profiling                     Profiling information from the search tree, including expansion models call and reactant generation
stock_info                    Dictionary of the stock availability for each of the starting material in all extracted routes
top_scores                    Comma-separated list of the score of the extracted routes (default to MCTS reward)
trees                         A list of the extracted routes as dictionaries
============================= ===========

If you running the tool with a single SMILES, all of this data will be printed to the screen, except
the ``stock_info`` and ``trees``.
