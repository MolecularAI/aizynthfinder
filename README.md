# aizynthfinder

aizynthfinder is a tool for retrosynthetic planning. The algorithm is based on a Monte Carlo tree search that recursively breaks down a molecule to purchasable precursors. The tree search is guided by a policy that suggests possible precursors by utilizing a neural network trained on a library of known reaction templates.  

## Prerequisites

Before you begin, ensure you have met the following requirements:

* Linux, Windows or Mac platforms are supported - as long as the dependencies are supported on these platforms.

The tool has been developed on a Linux platform, but the software has been tested on Windows 10 and macOS Catalina.

* You have installed anaconda or miniconda with python 3.6 or later


## Installation

To install aizynthfinder, follow these steps:

* First, install these conda packages

```
conda install -c rdkit rdkit -y
conda install -c anaconda tensorflow>=2.1.0 -y
conda install graphviz -y
```

if you have GPU and CUDA libraries enabled on your machine, you can install the ``tensorflow-gpu`` package instead.

* Secondly, install the ``aizynthfinder`` package

```
python -m pip install https://github.com/MolecularAI/aizynthfinder/archive/v1.0.0.tar.gz
```

if you want to install the latest version

or

```
python -m pip install -e .
```

if you are a developer, using the repository.

## Usage

The tool will install the ``aizynthcli`` and ``aizynthapp`` tools
as interfaces to the algorithm:

```
aizynthcli --config config_local.yml --smiles smiles.txt
aizynthapp --config config_local.yml
```

Consult the documentation here for more information.

To use tool you need

    1. A stock file
    2. A trained rollout policy network (including the Keras model and the list of unique templates)


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

* @SGenheden
* @EBjerrum
* @A-Thakkar

The contributors have limited time for support questions, but pleae do not hesitate to submit an issue (see above).

## License

The software is licensed under the MIT license (see LICENSE file), and is free and provided as-is.

## References

1. Thakkar A, Kogej T, Reymond J-L, et al (2019) Datasets and their influence on the development of computer assisted synthesis planning tools in the pharmaceutical domain. Chem Sci. https://doi.org/10.1039/C9SC04944D
