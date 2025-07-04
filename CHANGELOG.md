# CHANGELOG

## Version 4.4.0 2025-06-30

## Features

- Add scorers that are wrappers around rxnutils scorers
- Add filter policy for keeping substructures frozen
- Add stock option based on molecular weight

## Miscellaneous 

- Add rxnutils as dependency
- Remove routines in rxnutils to reduce code duplication

## Version 4.3.2 2025-04-07

### Trivial changes

- Upgraded dependencies

## Version 4.3.1 2024-12-05

### Trival changes

- Upgraded dependencies
- Removed conda env file for user setup
- Removed download of old HDF5 model files

### Bug-fixes

- Fixed import of clustering GUI causing issues

## Version 4.3.0 2024-05-22

### Features

- Focussed bonds feature has been created to allow users to choose bonds to freeze and break.
- Functionality has been added to run the broken bonds scorer with MCTS.
- A disconnection-aware Chemformer expansion strategy has been introduced in the plugins.
- Weights for single expansion policies can now be provided as input through the config file.
- The priors returned from the multi-expansion strategy have been rescaled.
- Functionality to mask reaction templates is now supported in the `TemplateBasedExpansionStrategy`.
- The `TemplateBasedDirectExpansionStrategy` has been implemented to directly apply the template in the search process.
- Added support for the C++ version of RDChiral.
- The multi-objective MCTS core algorithm is now implemented within the MCTS search functionality.
- A separate GUI component has been introduced for doing MO tree analysis. Additional functionalities have been added to GUI widgets to set two rewards/objectives to MCTS search. Pareto front is automatically plotted if MO-MCTS is run. Route re-ordering is automatically disabled for MO-MCTS.
- Preprocessing of the tree search can now be done using `aizynthcli`.
- The `StockAvailablityScorer` has been updated such that it takes and additional `other_source_score` parameter.
- A `cutoff_number` parameter can be provided to the multi-expansion strategy to obtain only the top predictions.
- A `BrokenBondsScorer` has been created for scoring nodes and reaction trees based on the breaking of atom bonds.
- A `RouteSimilarityScorer` has been created for scoring based on an LSTM model for computing Tree Edit Distance to a set of reference routes.
- A `DeltaSyntheticComplexityScorer` has been created for scoring nodes based on the delta-synthetic-complexity of the node and its parent 'horizon' steps up in the tree.

### Bug-fixes

- Fixed an issue of sending multiple fingerprints to the GPCR Tensorflow serving model.
- `aizynthcli` has been fixed after updating with multi-objective analysis such that the tool accurately informs the user if a target is solved or not.

## Version 4.0.0 2023-11-30

### Features

- `Configuration` now supports a `rescale_prior` property which rescales the priors in `TemplateBasedExpansionStrategy`.
- Functionality of `ScorerCollection` has been extended.
- Pricing is now supported in `InMemoryInchiKeyQuery`.
- Reward scorer has been added to Configuration scorers as `search_reward` item.
- `MaxTransformScorerer` and `FractionInStockScorer` have been created to separate scores in `StateScorer`.
- Reaction routes are scored with all reward scorers after the search is complete.
- `StockAvailabilityScorer` and `ReactionClassMembershipScorer` have been added to the scorers.
- Atom mapping existing in target molecule can be inherited.
- A caching feature has been add to the expansion strategies.
- Degenerate states can now be grouped in the MCTS algorithm.
- `ChemformerBasedExpansionStrategy` and `ModelZooExpansionStrategy` can be found under plugins and used as additional expansion strategies.
- Support for the stereocenter model has been added with to support chiral fingerprints for molecules and reactions.
- The `Configuration` format has been entirely revamped to a more easy-to-use format.
- A `MolbloomFilterQuery` has been created as a stock query class.

### Trivial changes

- Python version requirements have been updated to versions 3.9 - 3.11.
- `MoleculeCost` has been moved from aizynthfinder.context.cost.collection to the Retro* package.

### Deprecations

- Graphviz has been removed from aizynthfinder.utils.image.
- `Reaction` class has been removed from aizynthfinder.chem.reaction.
- aizynthfinder.context.cost package has been removed.

### Bug-fixes

- Fixed an issue with `max_transforms` to ensure only the given number of maximum depth is considered.
- Pinned Jupyter notebook version to ^6.5.3 to avoid errors when displaying widgets.
- Rollout child has been removed from the MCTS search algorithm.

## Version 3.7.0 2023-06-01 (2023-04-11)

### Features

- Environment variables can be loaded from YAML file.
- Restart feature in AiZynth CLI.
- Use ONNX models for expansion and filtering.
- Default filter is now loaded automatically in AiZynth CLI.

### Trivial changes

- Using isort to sort the imports of packages.
- Tensorflow is now an optional dependency.
- more-itertools dependencies removed.

### Breaking changes

- training modules has now been removed.

### Bug-fixes

- Fixing an issue when concatenating JSON data files.

## Version 3.6.0 2022-11-28 (2022-11-25)

### Features

- aizynthfinder can now be installed as a pure-python package
- ReactionTree now store at what iteration the route was created
- ReactionTree now supports a metadata property
- MCTS nodes populate a property that stores at what iteration the route was created
- aizynthcli now outputs stock information and route metadata and scores
- aizynthcli now has better error handling of invalid SMILES input
- Graphviz dependency for route drawing is removed, pure python implementation is used instead
- Dependencies reworked so that a minimal package can be installed
- Extra dependencies are related to route distances, clustering, training, MongoDB and external models
- Downloaded files are now the latest USPTO and Ringbreaker models

### Trivial changes

- RDKit and route-distances package now installed from pypi
- Documentation updated and extended

## Version 3.5.0 2022-11-28 (2022-07-21)

### Features

- Atom-mapping is tracked from product to reactant
- Support loading of template library from (gzipped) CSV file
- Support of saving aizynthcli output to (gzipped) JSON file
- AiZynthExpander now tracks non-applicable templates

### Bug-fixes
- Fixed failing test case

### Trival changes
- Silent progress bar when utilizing local Keras model

## Version 3.4.0 2022-04-28

### Features

- Gracefully fail predictions with aizynthcli if target is unsanitizable (Github issue 66)

### Trivial changes
- Update version of Sphinx, scikit-learn and Black
- Remove usage of depracted Pandas function append (Github issue 63)
- Correct documentation of aizynthcli (Github issue 67)

## Version 3.3.1 2022-03-09

### Trivial changes

- Updated pinned versions of route_distances, jinja2 and tensorflow

## Version 3.3.0 2022-02-24

### Features

- Support for Retro* tree search
- Support for breadth-first exhaustive tree search
- Support for depth-first proof-number tree search
- Possible to save concatenated reaction trees to separate file

### Bugfixes

- RouteCostScorer fix for rare routes

## Version 3.2.0 2022-02-24

### Features

- Profiling feature enabled in search trees
- New, customizable configuration of training pre-processing tools
- Generic post-processing support in aizynthcli
- Introduce short aliases for filter policies
- Reaction shape support in GraphViz visualisation of routes

## Version 3.1.0 2021-12-21

### Features

- ReactantsCountFiler (Github issue 42) to filter reactions with incompatible number of reactants
- ForwardRegenerationFiler to filter reactions where the forward reaction is incompatible
- Possible to skip quick Keras filter for specific policies
- Possible to select more than one policy in the GUI application
- Reaction classes has a hash function
- Possible to extract sub trees from ReactionTree objects
- RDKit can be used instead of RDChiral for expansions

### Bugfixes

- Possible to use more than depth 6 in the GUI application
- Fix failure in MctsNode class when expansion policy return no molecules

### Trivial changes

- Update type hints to be compatible with latest numpy release
- Update route-distances dependency

## Version 3.0.0  2021-07-26

### Features

- Generalized expansion and filter policies - policies no longer need to be Keras models
- Updated image generation for synthesis routes
- Improved code to prune "regeneration reactions" in MCTS (Github issue 38)
- Reactants are sorted in output from ReactionTree class
- New arguments to the "aizynthcli" tool
- Introduce option to return all solved routes from a TreeAnalysis class

### Breaking changes

- Package structure re-factorization
- Behaviour of RetroReaction, FilterPolicy and ExpansionPolicy classes has been changed
- Configuration class holding search tree settings is no longer loading settings from yaml-file on instantiation

### Deprecations

- Removed MCTS-specific routines in the TreeAnalysis class
- Removed code for identifying cycles in ReactionTree
- Remove JSON interface to AiZynthFinder class

### Bugfixes

- Fixed property assignment when converting MCTS node to ReactionTree

### Trivial changes

- Documentation updates
- Extensive re-factoring of test cases

## Version 2.6.0 - 2021-06-11 (2021-05-03)

### Features

- Add `AiZynthExpander` class as public interface to single-step reactions
- Route distance calculations and clustering is now dependent on package `route-distances`
- Route distance calculations with ML model is now supported
- Reaction tree objects now has property `is_branched`


## Version 2.5.0 - 2021-06-11 (2021-03-30)

### Features

- Introduce an interface for search algorithm working on AND/OR trees
- The progress bar is now removed upon completion of the tree search
- Move `ReactionTree` class to a new module (this will not brake backwards compatibility)

### Bug fixes

- Fix for Github issue #28 for a training script

### Trivial changes

- `scikit-learn` is now imported before `tensorflow`, according to Github issue 30

## Version 2.4.0 - 2021-02-22 (2021-02-22)


### Features

- Simplified interfaces for some internal classes
- Type hints have been added to the source code
- poetry is now used for dependency management

### Trivial changes

- Update of black version causing some re-formatting
- Small changes to documentation

## Version 2.3.0 - 2021-02-22 (2021-01-20)

### Features

- Add option to timeout distance calculation of route collections

## Version 2.2.1 - 2020-12-18

### Bugfixes

- Include templates for producing routes in package

## Version 2.2.0 - 2020-12-14 (2020-11-23)

### Features

- Add support for clustering of routes
- Add new scorers for total price and route score
- Add support for adding stop criteria to all stocks
- Add stop criteria UI for Jupyter notebook

### Trivial changes

- Refactor handling of MongoDB instances
- Refactor code to connect with REST and gRPC servers


## Version 2.1.0 - 2020-12-07 (2020-11-03)

### Features

- Possibility to extract a combined tree for all predicted routes
- New classes to use Tensorflow servers as policies

### Trivial changes

- Improved code for route visualizaition

## Version 2.0.0 - 2020-11-24 (2020-09-29)

### Breaking changes

- The `Reaction` class is renamed to `RetroReaction`
- The `Policy` class is renamed to `ExpansionPolicy`
- `config`, `stock`, `policy` and `scoring` modules moved to `aizynthfinder.context` package
- Tool "preprocess_rollout" is now called "preprocess_expansion"
- Some of the public methods of `Stock`, `Policy` and `ScorerCollection` classes are renamed
- Setting target molecule now destroys the search tree

### Features

- Add a filter policy to MCTS to remove unfeasible reactions
- Add tools to train filter policy
- Add logic to prevent cycle forming in MCTS by rejecting creation of parent molecule when expanding
- Introduce new `context` subpackage that contains the `config`, `stock`, `policy` and `scoring` modules
- The `Stock`, `ExpansionPolicy`, `FilterPolicy` and `ScorerCollection` classes now has a common interface for selection and loading
- Introduce possibility to remove unsantizable reactions from template library when training
- Catch exceptions from RDChiral more gracefully


## Version 1.2.0 - 2020-11-24 (2020-09-04)

### Features

- Add score "number of pre-cursors"
- Add score "number of pre-cursors in stock"
- Enable the loading of MCTS configuration from dict

### Bugfixes

- Fix bug in State score for branched reaction trees
- Fix bug in "Number of reaction" score for branched MCTS trees
- Fix bug in Jupyter GUI when selecting policies of several available

### Trivial changes

- Update python requirement specifications

## Version 1.1.0 - 2020-11-24 (2020-06-30)

### Features

- Introduce a route scoring framework
- Flag repetetive patterns in routes so they can be hidden
- Enable possibility to use more than one rollout policy
- Add --nproc argument to aizynthcli to enable trivial multiprocessing

### Bugfixes

- Add missing serialization of reaction metadata

### Trivial changes

- Add consistency check between rollout neural network and template list

## Version 1.0.0 - 2020-06-11

- First public version