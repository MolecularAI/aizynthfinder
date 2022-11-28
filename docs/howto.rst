How-to
=======

This page outlines a few guidelines on some more advanced use-cases of AiZynthFinder or
frequently raised issues.


Using Retro*
------------

AiZynthFinder implements other search algorithms than MCTS. This is an example of how Retro* can be used.

The search algorithm is specified in the configuration file.

.. code-block:: yaml

    properties:
      search_algorithm: aizynthfinder.search.retrostar.search_tree.SearchTree


This will use Retro* without a constant-valued oracle function. To specify the oracle function, you can 
do

.. code-block:: yaml

    properties:
      search_algorithm: aizynthfinder.search.retrostar.search_tree.SearchTree
    molecule_cost:
      aizynthfinder.search.retrostar.cost.RetroStarCost:
        model_path: retrostar_value_model.pickle

The pickle file can be downloaded from `here <https://github.com/MolecularAI/PaRoutes/blob/main/publication/retrostar_value_model.pickle?raw=true>`_ 


Using multiple expansion policies
---------------------------------

AiZynthFinder can use multiple expansion policies. This gives an example how a general USPTO and a RingBreaker model
can be used together

.. code-block:: yaml

    properties:
      additive_expansion: True
    policy:
      files:
        uspto:
          - uspto_keras_model.hdf5
          - uspto_unique_templates.csv.gz
        ringbreaker:
          - uspto_ringbreaker_keras_model.hdf5
          - uspto_ringbreaker_unique_templates.csv.gz


and then to use this with ``aizynthcli`` do something like this

.. code-block::

    aizynthcli --smiles smiles.txt --config config.yml --policy uspto ringbreaker


Output more routes
------------------

The number of routes in the output of ``aizynthcli`` can be controlled from the configuration file. 

This is how you can extract at least 25 routes but not more than 50 per target

.. code-block:: yaml

    properties:
        post_processing:
          min_routes: 25
          max_routes: 50

Alternatively, you can extract all solved routes. If a target is unsolved, it will return the number 
of routes specified by ``min_routes`` and ``max_routes``.

.. code-block:: yaml

    properties:
        post_processing:
          min_routes: 5
          max_routes: 10
          all_routes: True

