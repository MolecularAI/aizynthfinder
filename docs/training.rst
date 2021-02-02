Training
========

This page explains the training tools available in the `aizynthfinder` package.

Prerequisites
--------------

To start with, one needs a CSV file with templates, i.e. a pre-computed template
library. Tools to create such a library is not included in the `aizynthfinder` package,
because it is a tool that has a wider applicability than the tools provided in this package.
There are for instance tools like RdChiral (https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00286). 

The CSV file should contain the following columns in the given order:

* index - row index, this is ignored
* ID - an ID, this is ignored
* reaction_hash - a unique hash for each unique reaction. Are used to drop duplicates in the data.
* reactants - the SMILES of the reactants
* products - the SMILES of the products
* classification - reaction classification, this is used to add metadata to the tree search
* retro_template - the reaction template for the retro reaction
* template_hash - a unique hash for each unique template. Will be used to filter uncommon templates.
* selectivity - an indicator for selectivity, this is ignored
* outcomes - number of reaction outcomes, this is ignored

If you do not have all of these columns, or they are in another order, it is possible to
modify the ``library_headers`` setting in the configuration file (see below). However, the
``reaction_hash``, ``reactants``, ``products``, ``retro_template`` and ``template_hash`` 
columns are mandatory.

If you change this setting, you might also consider changing the ``metadata_headers`` that 
is a list of columns that is taken from the template library and injected into the tree search.
The default columns are "template_hash" and "classification", so if you for instance don't have
the ``classification`` column in your set, you need to update the ``metadata_headers``.

Training configuration
----------------------

The training is customizable to some extent by a configuration file in YAML format. If not provided,
the settings have default values. There are only a few settings that are of real interest to modify. 
They are shown in this snippet:

.. code-block:: yaml

    output_path: "."
    file_prefix: "full_uspto"
    batch_size: 256
    epochs: 100
    

These settings control the output directory, the prefix to all files,
the batch size for training and the number of training epochs.

For the other settings see the file ``default_training.yml`` in the ``data`` folder of the package. 

`Note!` The filename of the raw template library (discussed above) also needs to start with ``file_prefix`` and 
end of with ``_raw_template_library.csv``. 

Pre-processing and training
----------------------------

First the template library needs to be pre-processed such that

* The original library is pruned from templates that occur only a few times
* The template hash is turned into a label vector, i.e. the the target of the fitting
* The products SMILES are turned into fingerprint vectors, i.e. the input for the fitting
* The input and label matrices are split into training, testing and validation sets

This can be accomplished with


.. code-block:: 

    preprocess_expansion config.yaml


where ``config.yaml`` is your local configuration file for the training (see above).

Note that this procedure will take some time and might require a lot of memory.


Once this is done, you can train the network using

.. code-block::

    aizynth_training config.yaml expansion


Note that this might take a long time to converge.

The folder ``checkpoint`` will contain the Keras model that you can use as input 
to the tree search algorithm.

The pre-processing script created a file that ends with ``unique_templates.hdf5`` - 
this contains the unique templates and is the second input that you need for the tree search algorithm


Filter policy
-------------

To train a filter policy an array of tools are available. 

First, you need to generate negative data, i.e. reactions that are unfeasible. 

.. code-block::

    make_false_products config.yml strict
    make_false_products config.yml random
    make_false_products config.yml recommender


The first argument is a configuration file, similar to the one used above with the ``preprocess_expansion`` tool. 
The second argument should be "strict", "random" or "recommender" depending on what method you want to use.

When using the "recommender" method it is important to add the following to the configuration file:

.. code-block:: yaml

    recommender_model: "some_path/checkpoints/keras_model.hdf"

which points to the trained "recommender" model (see below).

The second step is pre-processing the training data:

.. code-block:: 

    preprocess_filter.py config.yaml


And the third and final step is the actual training:


.. code-block::

    aizynth_training config.yaml filter


The folder ``checkpoint`` will contain the Keras model that you can use as input 
to the tree search algorithm.


Training recommender model
--------------------------

Training to recommender model is very similar to training the expansion policy 


.. code-block:: 

    preprocess_recommender config.yaml
    aizynth_training config.yaml recommender


The folder ``checkpoint`` will contain the Keras model that you can use to generate negative data.