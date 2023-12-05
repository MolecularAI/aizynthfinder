
Configuration file
===================

The configuration file that is to be provided to the interfaces
specify file paths to policy models and stocks, in addition to
detailed parameters of the tree search.

Simple usage
------------

Let's say you have a

    * A trained Keras expansion model that is called `uspto_expansion.onnx`
    * A library of unique templates called `uspto_templates.csv.gz`
    * A stock file in HDF5 format, called `zinc_stock.hdf5`

(this could have been created by the ``smiles2stock`` tool, see :doc:`here <stocks>`)

then to use them in the tree search, create a file called `config.yml` that has the following content

.. code-block:: yaml

    expansion:
      full:
        - uspto_expansion.onnx
        - uspto_templates.csv.gz
    stock:
      zinc: zinc_stock.hdf5


Advanced usage
--------------

A more detailed configuration file is shown below

.. code-block:: yaml

    search:
      algorithm: mcts
      algorithm_config:
        C: 1.4
        default_prior: 0.5
        use_prior: True
        prune_cycles_in_search: True
        search_reward: state score
      max_transforms: 6
      iteration_limit: 100
      return_first: false
      time_limit: 120
      exclude_target_from_stock: True 
    expansion:
      my_policy:
        type: template-based
        model: /path/to/keras/model/weights.hdf5
        template: /path/to/hdf5/templates.hdf5
        template_column: retro_template
        cutoff_cumulative: 0.995
        cutoff_number: 50
        use_rdchiral: True
        use_remote_models: False
      my_full_policy:
        - /path/to/keras/model/weights.hdf5
        - /path/to/hdf5/templates.hdf5
    filter:
      uspto:
        type: quick-filter
        model: /path/to/keras/model/weights.hdf5
        exclude_from_policy: rc
        filter_cutoff: 0.05
        use_remote_models: False
      uspto_full: /path/to/keras/model/weights.hdf5
    stock:
      buyables:
        type: inchiset
        path: /path/to/stock1.hdf5
      emolecules: /path/to/stock1.hdf5

The (expansion) policy models are specified using two files
    * a checkpoint files from Keras in ONNX or hdf5 format,
    * a HDF5 or a CSV file containing templates.

A key like ``my_policy`` should be set and the configuration contains ``type``, ``model`` and ``template`` that must be provided. 
If the other settings are not assigned, their default values are taken. 
The policy ``my_full_policy`` exemplifies a short-cut to the template-based expansion model when no other settings need to be 
provided, only the ``model`` and ``templates`` can be provided. The default settings will be taken in this case.

The template file should be readable by ``pandas`` using  the ``table`` key and the ``retro_template`` column.
A policy can then be selected using the provided key, like ``my_policy`` in the above example.

The filter policy model is specified using a single checkpoint file.
Any key like ``uspto`` can be set. The settings contains ``type`` and ``model`` that must be provided. If the other
settings are not assigned, their default values are taken. 
The filter ``uspto_full`` exemplifies a short-cut to the quick-filter model when no other settings need to be 
provided, only the ``model`` can be provided. The default settings will be taken in this case.


The stock files can be
     * HDF5 files with the ``table`` key an the ``inchi_key`` column.
     * A CSV file with a ``inchi_key`` column
     * A text file a single column

In all cases, the column should contain pre-computed inchi keys of the molecules.
The stocks can be set using any key, like ``buyables`` or ``emolecules`` in the above example.
The ``type`` and ``path`` parameters can also be set along with other parameters.
If no other settings need to be provided, only the ``path`` can be provided, whereby it will be 
treated as a 
short-cut to the inchiset class..

The values in the ``search`` sections are optional, and if missing, default values are considered. These
values can also be taken from environment variables. An example of this can be seen as below:

.. code-block:: yaml

    search:
      iteration_limit: ${ITERATION_LIMIT}
      algorithm_config:
        C: ${C}
      time_limit: ${TIME_LIMIT}
      max_transforms: ${MAX_TRANSFORMS}

These are the available ``search`` settings. The ``algorithm_config`` refers to MCTS settings:

============================================ ============== ===========
Setting                                      Default value  Description
============================================ ============== ===========
algorithm                                    mcts           The search algorithm. Can be set to `package.module.ClassName` to use a custom search method.
algorithm_config: C                          1.4            The C value used to balance exploitation and exploration in the upper confidence bound score of the nodes.
algorithm_config: default_prior              0.5            The prior that is used if policy-provided priors are not used.
algorithm_config: use_prior                  True           If True, priors from the policy is used instead of the `default_prior`.
algorithm_config: prune_cycles_in_search     True           If True, prevents the MCTS from creating cycles by recreating previously seen molecules when it is expanded.
algorithm_config: search_reward              state score    The scoring used for the MCTS search algorithm.
algorithm_config: immediate_instantiation    []             list of expansion policies for which the MCTS algorithm immediately instantiate the children node upon expansion
algorithm_config: mcts_grouping              -              if is partial or full the MCTS algorithm will group expansions that produce the same state. If ``partial`` is used the equality will only be determined based on the expandable molecules, whereas ``full`` will check all molecules.
max_transforms                               6              The maximum depth of the search tree.
iteration_limit                              100            The maximum number of iterations for the tree search.
time_limit                                   120            The maximum number of seconds to complete the tree search.
return_first                                 False          If True, the tree search will be terminated as soon as one solution is found.
exclude_target_from_stock                    True           If True, the target is in stock will be broken down.
============================================ ============== ===========


The ``post_processing`` settings are:

============================================ ============== ===========
Setting                                      Default value  Description
============================================ ============== ===========
min_routes                                   5              The minumum number of routes to extract if ``all_routes`` is not set.
max_routes                                   25             The maximum number of routes to extract if ``all_routes`` is not set.
all_routes                                   False          If True, will extract all solved routes.
route_distance_model                         N/A            If set, will load the quick route distance model from this checkpoint file.
route_scorer                                 state score    The scoring for routes when extracting them in the post processing step.
============================================ ============== ===========


The ``expansion`` settings are for template-based models:

============================================ ============== ===========
Setting                                      Default value  Description
============================================ ============== ===========
template_column                              retro_template The column in the template file that contains the templates.
cutoff_cumulative                            0.995          The accumulative probability of the suggested templates is capped at this value. All other templates above this threshold are discarded.
cutoff_number                                50             The maximum number of templates that will be returned from the expansion policy.
use_rdchiral                                 True           If True, will apply templates with RDChiral, otherwise RDKit will be used.
use_remote_models                            False          If True, will try to connect to remote Tensorflow servers.
rescale_prior                                False          If True, will apply a softmax function to the priors.
============================================ ============== ===========


The ``filter`` settings are for quick-filter models:

============================================ ============== ===========
Setting                                      Default value  Description
============================================ ============== ===========
exclude_from_policy                          []             The list of names of the filter policies to exclude.
filter_cutoff                                0.05           The cut-off for the quick-filter policy.
use_remote_models                            False          If True, will try to connect to remote Tensorflow servers.
============================================ ============== ===========
