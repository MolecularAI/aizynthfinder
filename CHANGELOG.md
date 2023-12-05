# CHANGELOG

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

## Version 3.6.0 2022-11-25

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
- (AZ) support for template lookup is removed

### Trivial changes

- RDKit and route-distances package now installed from pypi
- Documentation updated and extended

## Version 3.5.2 2022-07-27

### Trivial changes

- Fix dependency for route-distances
- Update Seldon environment

## Version 3.5.1 2022-07-22

### Bug-fixes

- Fix atom-mapping inheritence for atoms only in reactants

## Version 3.5.0 2022-07-21

### Features

- All Seldon models now uses version 0.5.0 of aizynthfinder_models (Reaction connect models july 2022)
- Atom-mapping is tracked from product to reactant
- Support loading of template library from (gzipped) CSV file
- Support of saving aizynthcli output to (gzipped) JSON file
- DeepRetrosynthesis model now returns a field called "mapped_reaction_smiles"
- AiZynthExpander now tracks non-applicable templates

### Bug-fixes
- Fixed failing test case

### Trival changes
- Silent progress bar when utilizing local Keras model


## Version 3.4.0 2022-04-28

### Features

- Deep retrosynthesis Seldon model
- Gracefully fail predictions with aizynthcli if target is unsanitizable (Github issue 66)

### Trivial changes
- Update version of Sphinx, scikit-learn and Black
- Remove usage of depracted Pandas function append (Github issue 63)
- Correct documentation of aizynthcli (Github issue 67)

## Version 3.3.0 2022-02-14

### Features

- Support for breadth-first exhaustive tree search
- Support for depth-first proof-number tree search
- Possible to save concatenated reaction trees to separate file

### Trivial changes

- Retro* search algorithm moved to public package

## Version 3.2.1 2022-01-17

### Bugfixes

- RouteCostScorer fix for rare routes

## Version 3.2.0 2022-01-11

### Features

- Profiling feature enabled in search trees
- New, customizable configuration of training pre-processing tools
- Generic post-processing support in aizynthcli
- Introduce short aliases for filter policies
- Reaction shape support in GraphViz visualisation of routes

## Version 3.1.0 2021-11-09

### Features

- Reaction connect API is updated for version 2.0
- Filter to remove unreasonable expansions (when the template remove most of the molecule)
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

## Version 3.0.1 2021-08-31

### Trivial changes

- REST and GRPC Tensorflow clients now using named input vectors when provided 

### Bugfixes

- Correct input vectors provided to REST and GRPC tensorflow clients for filter policy

## Version 3.0.0  2021-07-26

### Features

- Generalized expansion and filter policies - policies no longer need to be Keras models
- ReactionConnect-based expansion policy and stock query
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

### Bugfixes

- Fixed property assignment when converting MCTS node to ReactionTree

### Trivial changes

- Documentation updates
- Extensive re-factoring of test cases

## Version 2.6.0 - 2021-05-03

### Features

- Add `AiZynthExpander` class as public interface to single-step reactions
- Route distance calculations and clustering is now dependent on package `route-distances`
- Route distance calculations with ML model is now supported
- Reaction tree objects now has property `is_branched`

### Trivial changes

- `scikit-learn` is now imported before `tensorflow`, according to Github issue 30

## Version 2.5.0 - 2021-03-30

### Features

- Introduce an interface for search algorithm working on AND/OR trees
- Add a reference implementation of the Retro* algorithm in the az sub-package
- Improved AiZynthExpander model to group predictions with the same precursors
- The progress bar is now removed upon completion of the tree search
- Move `ReactionTree` class to a new module (this will not brake backwards compatibility)

### Bug fixes

- Fix for Github issue #28 for a training script

### Trivial changes

- Use pylint as the code linter

## Version 2.4.0 - 2021-02-22

### Features

- Simplified interfaces for some internal classes
- Type hints have been added to the source code
- poetry is now used for dependency management

### Trivial changes

- Update of black version causing some re-formatting
- Small changes to documentation

## Version 2.3.0 - 2021-01-20

### Features

- Add option to timeout distance calculation of route collections
- Add option to perform clustering in Seldon AiZynthFinder model
- Add new Seldon model AiZynthExpander encapsulating expansion policies

### Bugfixes

- Jinja templates are now packaged with the python package

## Version 2.2.1 - 2020-12-01

### Bug fixes

- Fix price and amount function in AiZynth MongoDB stock query

### Trivial changes

- Introduce caching of fingerprints in template info class

## Version 2.2.0 - 2020-11-23

### Features

- Add support for clustering of routes
- Add new scorers for total price and route score
- Add support for adding stop criteria to all stocks
- Add stop criteria UI for Jupyter notebook
- Enable the applicability filtered-expansion model in Seldon model
- Improved memory management for template look-up

### Trivial changes

- Refactor handling of MongoDB instances
- Refactor code to connect with REST and gRPC servers
- Refactor AiZynth buyables query class


## Version 2.1.1 - 2020-11-08

### Bugfixes

- Matching names of policies and scorers in Seldon model

## Version 2.1.0 - 2020-11-03

### Features

- Possibility to extract a combined tree for all predicted routes
- New classes to use Tensorflow servers as policies
- Add an internal stock query class to use on the AiZynth platform
- Enable use of the tool as a Seldon server

### Trivial changes

- Reduced memory consumption for average product similarity scorer
- Improved code for route visualizaition

## Version 2.0.0 - 2020-09-29

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


## Version 1.2.0 - 2020-09-04

### Features

- Add score "number of pre-cursors"
- Add score "number of pre-cursors in stock"
- Add a stock query to use the AiZynth Buyable API
- Enable the loading of MCTS configuration from dict

### Bugfixes

- Fix bug in State score for branched reaction trees
- Fix bug in "Number of reaction" score for branched MCTS trees
- Fix bug in "Average template occurence" score for trees with no reactions
- Fix bug in Jupyter GUI when selecting policies of several available

### Trivial changes

- Update python requirement specifications

## Version 1.1.0 - 2020-06-30

### Features

- Introduce a route scoring framework
- Flag repetetive patterns in routes so they can be hidden
- Enable possibility to use more than one rollout policy
- Add --nproc argument to aizynthcli to enable trivial multiprocessing
- Add tool to download publicly released data

### Bugfixes

- Add missing serialization of reaction metadata

### Trivial changes

- Add consistency check between rollout neural network and template list

## Version 1.0.0 - 2020-06-01

### Breaking changes

- The tree search algorithm now breaks down any target molecule in stock by default

### Features

- Add option to break down target molecules that are in stock
- Add a GUI extension and API to lookup reactions in Reaxys and Pistachio
- Metadata is added to reactions and to the output

### Bugfixes

- Fix issue with tree search when all children are rejected 

### Trivial changes

- Improved documentation
- Improved test coverage
- Interface to graphviz is entirely within the package, not with third-party library

## Version 0.3.0 - 2020-04-14

### Features

- Add tools to train rollout policy model
- Add tools to create stock from SMILES
- Introduced the RouteCollection class for easy analysis of routes
- Introduced possibility to use external Keras model

### Bugfixes

- Fix recreation of reaction tree from dict for certain types of routes


## Version 0.2.0 - 2020-03-16

### Features

- Moved internal tree search classes to aizynthfinder.mcts package
- Removed concept of runaway nodes, now exhaustively look through all actions
- Add auto-selecting of stocks and policy in aizynthcli tool
- Add support for custom stock class in aizynthcli tool
- Add support for MongoDB stock lookup

### Bugfixes

- Fix cut-off of predictions in Policy class based on threshold and minimum number

## Version 0.1.0 - 2020-02-18

- Re-factored codebase
- Provide route visualization in Jupyter Notebooks
- Provide command to run predictions in batch