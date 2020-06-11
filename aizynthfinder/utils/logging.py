""" Module containing routines to setup proper logging
"""
import logging.config
import os
import sys

import yaml
import tensorflow
from rdkit import RDLogger

from aizynthfinder.utils.paths import data_path


# This supress the printing of the Keras backend
stderr = sys.stderr
sys.stderr = open(os.devnull, "w")
import tensorflow.keras  # noqa

sys.stderr = stderr

# Supress tensforflow logging
tf_logger = tensorflow.get_logger()
tf_logger.setLevel(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

# Suppress RDKit errors due to incomplete template (e.g. aromatic non-ring atoms)
rd_logger = RDLogger.logger()
rd_logger.setLevel(RDLogger.CRITICAL)


def logger():
    """
    Returns the logger that should be used by all classes

    :return: the logger object
    :rtype: logging.Logger
    """
    return logging.getLogger("aizynthfinder")


def setup_logger(console_level, file_level=None):
    """
    Setup the logger that should be used by all classes

    The logger configuration is read from the `logging.yml` file.

    :param console_level: the level of logging to the console
    :type console_level: int
    :param file_level: the level of logging to file, if not set logging to file is disabled, default to None
    :type file_level: int, optional
    :return: the logger object
    :rtype: logging.Logger
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
