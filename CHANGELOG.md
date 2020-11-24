# CHANGELOG

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