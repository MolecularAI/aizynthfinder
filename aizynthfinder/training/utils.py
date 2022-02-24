""" Module containing various classes and routines used in training tools
"""
from __future__ import annotations
import os
import hashlib
from collections.abc import Mapping
from typing import TYPE_CHECKING

import yaml
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import train_test_split


from aizynthfinder.utils.paths import data_path
from aizynthfinder.chem import Molecule, MoleculeException

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        Any,
        List,
        Union,
        Tuple,
        Sequence,
    )


class Config:
    """
    Class that holds the configuration of the training.

    The settings are read from a yaml file. Default values
    for all settings are taken from the ``data`` directory of the repo.

    Settings can be read and set with

    .. code-block::

        batch_size = config["batch_size"]
        config["batch_size"] = 100

    :param config_filename: the path to a yaml file with settings
    """

    def __init__(self, config_filename: str = None) -> None:
        filename = os.path.join(data_path(), "default_training.yml")
        with open(filename, "r") as fileobj:
            default_config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)

        self._config = default_config

        if config_filename is None:
            return

        with open(config_filename, "r") as fileobj:
            user_config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
        self._update_dict(default_config, user_config)

    def __getitem__(self, item: str) -> Any:
        return self._config[item]

    def __setitem__(self, item: str, value: Any) -> None:
        self._config[item] = value

    def filename(self, label: str) -> str:
        """
        Return the absolute path to a file specified partly
        by settings.

        A filename is constructed from:
        ``output_path + file_prefix + file_postfix``

        where ``file_postfix`` is either taken from the settings
        by look-up, e.g. ``config["file_postfix"][label]`` or
        directly from the ``label`` argument.

        :param label: the file postfix
        :return: the filepath
        """
        return os.path.join(
            self["output_path"],
            self["file_prefix"] + self["file_postfix"].get(label, label),
        )

    @staticmethod
    def _update_dict(original: StrDict, other: StrDict) -> StrDict:
        # Used to complement the update method of the built-in dict type
        # it works for recursive dicts (to 1 level)
        for key, val in original.items():
            if key not in other or not isinstance(other[key], type(val)):
                continue
            if isinstance(val, Mapping):
                original[key] = Config._update_dict(original[key], other[key])
            else:
                original[key] = other[key]
        for key, val in other.items():
            if key not in original:
                original[key] = val
        return original


def create_reactants_molecules(reactants_str: str) -> List[Molecule]:
    """
    Create Molecule objects from a SMILE string of reactants.

    Only molecules with atom mapping is kept.

    :param reactants_str: the SMILES string of the reactants
    :return: the Molecule objects
    """
    mols = []
    for smiles in reactants_str.split("."):
        try:
            mol = Molecule(smiles=smiles, sanitize=True)
        except MoleculeException:
            pass
        else:
            if mol.has_atom_mapping():
                mols.append(mol)
    return mols


def is_sanitizable(args: Tuple[str]) -> bool:
    """
    Check whether a SMILES is sanitizable
    :param args: the SMILES in the first element
    :return: whether the SMILES is sanitizable
    """
    smiles = args[0]
    try:
        Molecule(smiles=smiles, sanitize=True)
    except MoleculeException:
        return False
    else:
        return True


def reverse_template(retro_template: str) -> str:
    """
    Reverse the reaction template to switch product and reactants

    :param retro_template: the reaction template
    :return: the reverse template
    """
    return ">>".join(retro_template.split(">>")[::-1])


def reaction_hash(reactants_smiles: str, product: Molecule) -> str:
    """
    Create a reaction hash

    :param reactants_smiles: the SMILES string of the reactants
    :param product: the product molecule
    :return: the hash
    """
    reactant_inchi = Molecule(smiles=reactants_smiles).inchi
    product_inchi = product.inchi
    concat_inchi = reactant_inchi + "++" + product_inchi
    return hashlib.sha224(concat_inchi.encode("utf8")).hexdigest()


def split_and_save_data(
    data: Union[pd.DataFrame, np.ndarray, sparse.csr_matrix],
    data_label: str,
    config: Config,
) -> None:
    """
    Split input data into training, testing and validation sets,
    and then saves it to disc.

    The input data can be either a pandas DataFrame a numpy array
    or a sparse matrix.

    :param data: the data to split
    :param data_label: the label of the data, if its input or labels
    :param config: the settings
    """
    train_size = config["split_size"]["training"]
    testing_frac = config["split_size"]["testing"]
    validation_frac = config["split_size"]["validation"]
    testing_size = testing_frac / (testing_frac + validation_frac)

    train_arr, test_arr = train_test_split(
        data, train_size=train_size, random_state=42, shuffle=True
    )
    val_arr, test_arr = train_test_split(
        test_arr, test_size=testing_size, random_state=42, shuffle=True
    )

    array_dict = {"training_": train_arr, "validation_": val_arr, "testing_": test_arr}
    for label_prefix, arr in array_dict.items():
        filename = config.filename(label_prefix + data_label)
        if isinstance(data, pd.DataFrame):
            arr.to_csv(
                filename,
                sep=config["csv_sep"],
                header=config["in_csv_headers"],
                index=False,
            )
        elif isinstance(data, np.ndarray):
            np.savez(filename, arr)
        else:
            sparse.save_npz(filename, arr, compressed=True)


def smiles_to_fingerprint(args: Sequence[str], config: Config) -> np.ndarray:
    """
    Convert a SMILES to a fingerprint vector

    :param args: the SMILES in the first element
    :param config: the settings
    :return: the fingerprint
    """
    smiles = args[0]
    return (
        Molecule(smiles=smiles)
        .fingerprint(
            config["fingerprint_radius"],
            config["fingerprint_len"],
        )
        .astype(np.int8)
    )


def reactants_to_fingerprint(args: Sequence[str], config: Config) -> np.ndarray:
    """
    Convert a SMILES string of reactants to a fingerprint

    :param args: the SMILES in the first element
    :param config: the settings
    :return: the fingerprint
    """
    reactants_smiles = args[0]
    fingerprints = []
    for smiles in reactants_smiles.split("."):
        try:
            mol = Molecule(smiles=smiles, sanitize=True)
        except MoleculeException:
            pass
        else:
            if mol.has_atom_mapping():
                fingerprints.append(
                    mol.fingerprint(
                        config["fingerprint_radius"], config["fingerprint_len"]
                    )
                )
    return sum(fingerprints)  # type: ignore


def reaction_to_fingerprints(args: Sequence[str], config: Config) -> np.ndarray:
    """
    Convert a reaction SMILEs string  a fingerprint

    :param args: the product SMILES in the first element, and reactants SMILES in the second
    :param config: the settings
    :return: the fingerprint
    """
    product_smiles, reactants_smiles = args
    product_fp = smiles_to_fingerprint([product_smiles], config)
    reactant_fp = reactants_to_fingerprint([reactants_smiles], config)

    return (product_fp - reactant_fp).astype(np.int8)


def split_reaction_smiles(data: pd.DataFrame, config: Config) -> pd.DataFrame:
    """
    Split a column of reaction SMILES into reactant and product columns

    :param data: the dateframe to process
    :param config: the training configuration
    :return: the new dataframe
    """
    smiles_split = data[config["reaction_smiles_column"]].str.split(">", expand=True)
    return data.assign(
        **{
            config["column_map"]["reactants"]: smiles_split[0],
            config["column_map"]["products"]: smiles_split[2],
        }
    )
