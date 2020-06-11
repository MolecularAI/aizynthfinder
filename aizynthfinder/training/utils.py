""" Module containing various classes and routines used in training tools
"""
import os
from collections.abc import Mapping

import yaml
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import train_test_split


from aizynthfinder.utils.paths import data_path
from aizynthfinder.chem import Molecule


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
    :type config_filename: str, optional
    """

    def __init__(self, config_filename=None):
        filename = os.path.join(data_path(), "default_training.yml")
        with open(filename, "r") as fileobj:
            default_config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)

        self._config = default_config

        if config_filename is None:
            return

        with open(config_filename, "r") as fileobj:
            user_config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
        self._update_dict(default_config, user_config)

    def __getitem__(self, item):
        return self._config[item]

    def __setitem__(self, item, value):
        self._config[item] = value

    def filename(self, label):
        """
        Return the absolute path to a file specified partly
        by settings.

        A filename is constructed from:
        ``output_path + file_prefix + file_postfix``

        where ``file_postfix`` is either taken from the settings
        by look-up, e.g. ``config["file_postfix"][label]`` or
        directly from the ``label`` argument.

        :param label: the file postfix
        :type label: str
        :return: the filepath
        :rtype: str
        """
        return os.path.join(
            self["output_path"],
            self["file_prefix"] + self["file_postfix"].get(label, label),
        )

    @staticmethod
    def _update_dict(original, other):
        for key, val in original.items():
            if key not in other or not isinstance(other[key], type(val)):
                continue
            if isinstance(val, Mapping):
                original[key] = Config._update_dict(original[key], other[key])
            else:
                original[key] = other[key]
        return original


def split_and_save_data(data, data_label, config):
    """
    Split input data into training, testing and validation sets,
    and then saves it to disc.

    The input data can be either a pandas DataFrame a numpy array
    or a sparse matrix.

    :param data: the data to split
    :type data: pandas.DataFrame or np.ndarray or scipy.sparse object
    :param data_label: the label of the data, if its input or labels
    :type data_label: str
    :param config: the settings
    :type config: Config
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
            arr.to_csv(filename, mode="w", header=False, index=False)
        elif isinstance(data, np.ndarray):
            np.savez(filename, arr)
        else:
            sparse.save_npz(filename, arr, compressed=True)


def smiles_to_fingerprint(args, config):
    """
    Convert a SMILES to a fingerprint vector

    :param args: the SMILES in the first element
    :type args: tuple
    :param config: the settings
    :type config: Config
    :return: the fingerprint
    :rtype: numpy.ndarray
    """
    smiles = args[0]
    return Molecule(smiles=smiles).fingerprint(
        config["fingerprint_radius"], config["fingerprint_len"],
    )
