
Configuration file
===================

The configuration file that is to be provided to the interfaces 
specify file paths to policy models and stocks, in addition to 
detailed parameters of the tree search. 

Simple usage
------------

Let say you have a 

    * A trained Keras expansion model that is called `full_uspto_model.hdf5`
    * A library of unique templates called `full_uspto_templates.hdf5`

(these could have been created by the training tools, see :doc:`here <training>`)

    * A stock file in HDF5 format, called `zinc_stock.hdf5`

(this could have been created by the ``smiles2stock`` tool, see :doc:`here <stocks>`)

then to use them in the tree search, create a file called `config.yml` that has the following content 

.. code-block:: yaml

    policy:
      files:
        full_uspto:
          - full_uspto_model.hdf5
          - full_uspto_templates.hdf5
    stock:
      files:
        zinc: zinc_stock.hdf5


Advanced usage
--------------

A more detailed configuration file is shown below

.. code-block:: yaml

    properties:
      iteration_limit: 100
      return_first: false
      time_limit: 120
      C: 1.4
      cutoff_cumulative: 0.995
      cutoff_number: 50
      max_transforms: 6
    policy:
      files:
        my_policy:
          - /path/to/keras/model/weights.hdf5
          - /path/to/hdf5/templates.hdf5
    filter:
      files:
        my_policy: /path/to/keras/model/weights.hdf5
    stock:
      files:
        stock1: /path/to/stock1.hdf5
        stock2: /path/to/stock1.hdf5

The (expansion) policy models are specified using two files
    * a checkpoint files from Keras in hdf5 format,
    * a HDF5 file containing templates. 

The filter policy model is specified using a single checkpoint file from Keras in hdf5 format.
    
The template file should be readable by ``pandas`` using  the ``table`` key and the ``retro_template`` column. 
A policy can then be selected using the provided key, like ``my_policy`` in the above example. 

The stock files should be HDF5 files with the ``table`` key an the ``inchi_key`` column. The column
should contain pre-computed inchi keys of the molecules. 
The stocks can then be selected using the provided key, like ``stock1`` or ``stock2`` in the above example.

The values in the ``properties`` sections are optional, and if missing, default values are provided.

These are the available properties:

========================= ============== ===========
Property                  Default value  Description
========================= ============== ===========
C                         1.4            The C value used to balance exploitation and exploration in the upper confidence bound score of the nodes.
cutoff_cumulative         0.995          The accumulative probability of the suggested templates is capped at this value. All other templates above this threshold are discarded. 
cutoff_number             50             The maximum number of templates that will be returned from the expansion policy.
max_transforms            6              The maximum depth of the search tree
default_prior             0.5            The prior that is used if policy-provided priors are not used
use_prior                 True           If true, priors from the policy is used instead of the `default_prior`
return_first              False          If true, the tree search will be terminated as soon as one solution is found
iteration_limit           100            The maximum number of iterations for the tree search
time_limit                120            The maximum number of seconds to complete the tree search
exclude_target_from_stock True           If the target is in stock it will be broken down if this property is True
template_column           retro_template the column in the template file that contains the templates
filter_cutoff             0.05           the cut-off for the quick-filter policy
prune_cycles_in_search    True           prevents the MCTS from creating cycles by recreating previously seen molecules when it is expanded
additive_expansion        False          If true, reactions from all selected expansion policies will be appended, otherwise only the first non-empty expansion will be used
search_algorithm          mcts           The search algorithm. Can be set to `package.module.ClassName` to use a custom search method
========================= ============== ===========