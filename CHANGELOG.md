# CHANGELOG

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