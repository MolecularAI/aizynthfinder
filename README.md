# aizynthfinder

aizynthfinder is a tool for retrosynthetic planning. The algorithm is based on a Monte Carlo tree search that recursively breaks down a molecule to purchasable precursors. The tree search is guided by a policy that suggests possible precursors by utilizing a neural network trained on a library of known reaction templates.  

## Prerequisites

Before you begin, ensure you have met the following requirements:

* Linux, Windows or Mac platforms are supported - as long as the dependencies are supported on these platforms.

The tool has been developed on a Linux platform, but the software has been tested on Windows 10 and macOS Catalina.

* You have installed anaconda or miniconda with python 3.6, 3.7 or 3.8


## Installation

To install aizynthfinder, follow these steps:

* First, install these conda packages

```
conda install -c conda-forge "rdkit=>2019.09.1" -y
conda install -c anaconda "tensorflow>=2.1.0" -y
conda install graphviz -y
```

* Secondly, install the ``aizynthfinder`` package

```
python -m pip install https://github.com/MolecularAI/aizynthfinder/archive/v2.2.1.tar.gz
```

if you want to install the latest version

or

```
python -m pip install -e .
```

if you are a developer, using the repository.

Note on the graphviz installation: this package does not depend on any third-party python interfaces to graphviz but instead calls the `neato` and `dot` executables directly. If these executable are not in the `$PATH` environmental variable, the generation of route images will not work. If unable to install it properly with the default conda chanel, try using `-c anaconda`. 

## Usage

The tool will install the ``aizynthcli`` and ``aizynthapp`` tools
as interfaces to the algorithm:

```
aizynthcli --config config.yml --smiles smiles.txt
aizynthapp --config config.yml
```

Consult the documentation [here](https://molecularai.github.io/aizynthfinder/) for more information.

To use the tool you need

    1. A stock file
    2. A trained rollout policy network (including the Keras model and the list of unique templates)
    3. A trained filer policy network (optional)

Such files can be downloaded from [figshare](https://figshare.com/articles/AiZynthFinder_a_fast_robust_and_flexible_open-source_software_for_retrosynthetic_planning/12334577) and [here](https://figshare.com/articles/dataset/A_quick_policy_to_filter_reactions_based_on_feasibility_in_AI-guided_retrosynthetic_planning/13280507) or they can be downloaded automatically using

```
download_public_data my_folder
```

where ``my_folder`` is the folder that you want download to. 
This will create a ``config.yml`` file that you can use with either ``aizynthcli`` or ``aizynthapp``.

### Testing

Tests uses the ``pytest`` package. 

To use, first install the dependencies

```
python -m pip install -r requirements_dev.txt
```

and then run the tests using

```
pytest -v
```

## Contributing

We welcome contributions, in the form of issues or pull requests.

If you have a question or want to report a bug, please submit an issue.


To contribute with code to the project, follow these steps:

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Make your changes and commit them: `git commit -m '<commit_message>'`
4. Push to the remote branch: `git push`
5. Create the pull request.

Please use ``black`` package for formatting, and follow ``pep8`` style guide.


## Contributors

* [@SGenheden](https://www.github.com/SGenheden)
* [@EBjerrum](https://www.github.com/EBjerrum)
* [@A-Thakkar](https://www.github.com/A-Thakkar)

The contributors have limited time for support questions, but please do not hesitate to submit an issue (see above).

## License

The software is licensed under the MIT license (see LICENSE file), and is free and provided as-is.

## References

1. Thakkar A, Kogej T, Reymond J-L, et al (2019) Datasets and their influence on the development of computer assisted synthesis planning tools in the pharmaceutical domain. Chem Sci. https://doi.org/10.1039/C9SC04944D
2. Genheden S, Thakkar A, Chadimova V, et al (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. J. Cheminf. https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00472-1
