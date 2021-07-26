""" Module containing routines to setup proper logging
"""
# pylint: disable=ungrouped-imports, wrong-import-order, wrong-import-position, unused-import
import logging.config
import os

import yaml

# See Github issue 30 why sklearn is imported here
import sklearn  # noqa
from rdkit import RDLogger

from aizynthfinder.utils.paths import data_path

# Suppress tensforflow logging
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

import tensorflow  # noqa

tf_logger = tensorflow.get_logger()
tf_logger.setLevel(logging.WARNING)

# Suppress RDKit errors due to incomplete template (e.g. aromatic non-ring atoms)
rd_logger = RDLogger.logger()
rd_logger.setLevel(RDLogger.CRITICAL)


def logger() -> logging.Logger:
    """
    Returns the logger that should be used by all classes

    :return: the logger object
    """
    return logging.getLogger("aizynthfinder")


def setup_logger(console_level: int, file_level: int = None) -> logging.Logger:
    """
    Setup the logger that should be used by all classes

    The logger configuration is read from the `logging.yml` file.

    :param console_level: the level of logging to the console
    :param file_level: the level of logging to file, if not set logging to file is disabled, default to None
    :return: the logger object
    """
    filename = os.path.join(data_path(), "logging.yml")
    with open(filename, "r") as fileobj:
        config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)

    config["handlers"]["console"]["level"] = console_level
    if file_level:
        config["handlers"]["file"]["level"] = file_level
    else:
        del config["handlers"]["file"]
        config["loggers"]["aizynthfinder"]["handlers"].remove("file")

    logging.config.dictConfig(config)
    return logger()
