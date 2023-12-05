# Plugins

This folder contains some features for `aizynthfinder` that
does not yet fit into the main codebase. It could be experimental
features, or features that require the user to install some
additional third-party dependencies.

For the expansion models, you generally need to add the `plugins` folder to the `PYTHONPATH`, e.g.

    export PYTHONPATH=~/aizynthfinder/plugins/

where the `aizynthfinder` repository is in the home folder

## Chemformer expansion model

An expansion model using a REST API for the Chemformer model
is supplied in the `expansion_strategies` module.

To use it, you first need to install the `chemformer` package
and launch the REST API service that comes with it.

To use the expansion model in `aizynthfinder` you can use a config-file
containing these lines

    expansion:
        chemformer:
            type: expansion_strategies.ChemformerBasedExpansionStrategy
            url: http://localhost:8000/chemformer-api/predict
    search:
        algorithm_config:
            immediate_instantiation: [chemformer]
        time_limit: 300

The `time_limit` is a recommandation for allowing the more expensive expansion model
to finish a sufficient number of retrosynthesis iterations.

You would have to change `localhost:8000` to the name and port of the machine hosting the REST service.

You can then use the config-file with either `aizynthcli` or the Jupyter notebook interface.

## ModelZoo expansion model

An expansion model using the ModelZoo feature is supplied in the `expansion_strategies` 
module. This is an adoption of the code from this repo: `https://github.com/AlanHassen/modelsmatter` that were used in the publications [Models Matter: The Impact of Single-Step Models on Synthesis Prediction](https://arxiv.org/abs/2308.05522) and [Mind the Retrosynthesis Gap: Bridging the divide between Single-step and Multi-step Retrosynthesis Prediction](https://openreview.net/forum?id=LjdtY0hM7tf).

To use it, you first need to install the `modelsmatter_modelzoo` package from
https://github.com/PTorrenPeraire/modelsmatter_modelzoo and set up the `ssbenchmark` 
environment. 

Ensure that the `external_models` sub-package contains the models required.
If it does not, you will need to manually clone the required model repositories 
within `external_models`.

To use the expansion model in `aizynthfinder`, you can specify it in the config-file
under `expansion`. Here is an example setting to use the expansion model with `chemformer` 
as the external model:

    expansion:
        chemformer:
            type: expansion_strategies.ModelZooExpansionStrategy:
            module_path: /path_to_folder_containing_cloned_repository/modelsmatter_modelzoo/external_models/modelsmatter_chemformer_hpc/
            use_gpu: False
            params:
                module_path: /path_to_model_file/chemformer_backward.ckpt
                vocab_path: /path_to_vocab_file/bart_vocab_downstream.txt
    search:
        algorithm_config:
            immediate_instantiation: [chemformer]
        time_limit: 300