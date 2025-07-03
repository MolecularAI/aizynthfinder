How-to
=======

This page outlines a few guidelines on some more advanced use-cases of AiZynthFinder or
frequently raised issues.


Using Retro*
------------

AiZynthFinder implements other search algorithms than MCTS. This is an example of how Retro* can be used.

The search algorithm is specified in the configuration file.

.. code-block:: yaml

    search:
      algorithm: aizynthfinder.search.retrostar.search_tree.SearchTree


This will use Retro* without a constant-valued oracle function. To specify the oracle function, you can 
do

.. code-block:: yaml

    search:
      algorithm: aizynthfinder.search.retrostar.search_tree.SearchTree
      algorithm_config:
        molecule_cost: 
          cost: aizynthfinder.search.retrostar.cost.RetroStarCost
          model_path: retrostar_value_model.pickle
          fingerprint_length: 2048
          fingerprint_radius: 2
          dropout_rate: 0.1

The pickle file can be downloaded from `here <https://github.com/MolecularAI/PaRoutes/blob/main/publication/retrostar_value_model.pickle?raw=true>`_ 


Using multiple expansion policies
---------------------------------

AiZynthFinder can use multiple expansion policies. This gives an example how a general USPTO and a RingBreaker model
can be used together

.. code-block:: yaml

    expansion:
      uspto:
        - uspto_keras_model.hdf5
        - uspto_unique_templates.csv.gz
      ringbreaker:
        - uspto_ringbreaker_keras_model.hdf5
        - uspto_ringbreaker_unique_templates.csv.gz
      multi_expansion_strategy:
        type: aizynthfinder.context.policy.MultiExpansionStrategy
        expansion_strategies: [uspto, ringbreaker]
        additive_expansion: True

and then to use this with ``aizynthcli`` do something like this

.. code-block::

    aizynthcli --smiles smiles.txt --config config.yml --policy multi_expansion_strategy


Output more routes
------------------

The number of routes in the output of ``aizynthcli`` can be controlled from the configuration file. 

This is how you can extract at least 25 routes but not more than 50 per target

.. code-block:: yaml

    post_processing:
      min_routes: 25
      max_routes: 50

Alternatively, you can extract all solved routes. If a target is unsolved, it will return the number 
of routes specified by ``min_routes`` and ``max_routes``.

.. code-block:: yaml

    post_processing:
        min_routes: 5
        max_routes: 10
        all_routes: True


Running multi-objective (MO) MCTS with disconnection-aware Chemformer
------------------
Disconnection-aware retrosynthesis can be done with 1) MO-MCTS (state score + broken bonds score), 2) Chemformer or 3) both.

First, you need to specify the bond constraints under search, see below.
To run the MO-MCTS with the "broken bonds" score, add the "broken bonds" score to the list of search_rewards:

.. code-block:: yaml

  search:
      break_bonds: [[1, 2], [3, 4]]
      freeze_bonds: []
      algorithm_config:
        search_rewards: ["state score", "broken bonds"]

To use the disconnection-aware Chemformer, you first need to add the `plugins` folder to the `PYTHONPATH`, e.g.

    export PYTHONPATH=~/aizynthfinder/plugins/

The script for starting a disconnection-aware Chemformer service is available at https://github.com/MolecularAI/Chemformer.
The multi-expansion policy with template-based model and Chemformer is specified with:

.. code-block:: yaml
  expansion:
    standard:
      type: template-based
      model: path/to/model
      template: path/to/templates
    chemformer_disconnect:
      type: expansion_strategies.DisconnectionAwareExpansionStrategy
      url: "http://localhost:8023/chemformer-disconnect-api/predict-disconnection"
      n_beams: 5
    multi_expansion:
      type: aizynthfinder.context.policy.MultiExpansionStrategy
      expansion_strategies: [chemformer_disconnect, standard]
      additive_expansion: True
      cutoff_number: 50


To use MO-tree ranking and building, set:

.. code-block:: yaml
  post_processing:
    route_scorers: ["state score", "broken bonds"]
    
Note: If post_processing.route_scorers is not specified, it will default to search.algorithm_config.search_rewards.