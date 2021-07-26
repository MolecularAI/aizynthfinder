# CHANGELOG

## Version 3.0.0  2020-07-26

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

## Version 2.6.0 - 2020-06-11 (2020-05-03)

### Features

- Add `AiZynthExpander` class as public interface to single-step reactions
- Route distance calculations and clustering is now dependent on package `route-distances`
- Route distance calculations with ML model is now supported
- Reaction tree objects now has property `is_branched`


## Version 2.5.0 - 2020-06-11 (2020-03-30)

### Features

- Introduce an interface for search algorithm working on AND/OR trees
- The progress bar is now removed upon completion of the tree search
- Move `ReactionTree` class to a new module (this will not brake backwards compatibility)

### Bug fixes

- Fix for Github issue #28 for a training script

### Trivial changes

- `scikit-learn` is now imported before `tensorflow`, according to Github issue 30

## Version 2.4.0 - 2020-02-22 (2020-02-22)


### Features

- Simplified interfaces for some internal classes
- Type hints have been added to the source code
- poetry is now used for dependency management

### Trivial changes

- Update of black version causing some re-formatting
- Small changes to documentation

## Version 2.3.0 - 2020-02-22 (2020-01-20)

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