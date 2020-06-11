""" Module containing a class for encapsulating the settings of the tree search
"""
import os

import yaml

from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.paths import data_path
from aizynthfinder.mcts.policy import Policy
from aizynthfinder.mcts.stock import Stock, MongoDbInchiKeyQuery


class Configuration:
    """
    Encapsulating the settings of the tree search, including the policy,
    the stock and various parameters.

    All the parameters can be retrieved as attributes of the Configuration
    object, e.g.

    .. code-block::

        config.max_transforms  # The maximum of transform
        config.iteration_limit  # The maximum number of iterations


    On instantiation it will read default parameters from a config.yml
    file located in the `data` folder of the package.
    """

    def __init__(self):
        self._properties = {}
        filename = os.path.join(data_path(), "config.yml")
        with open(filename, "r") as fileobj:
            _config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
        self._update_from_config(_config)

        self.stock = Stock()
        self.policy = Policy(self)
        self._logger = logger()

    def __eq__(self, other):
        return self._properties == other._properties

    @classmethod
    def from_file(cls, filename):
        """
        Loads a configuration from a yaml file.
        The parameters not set in the yaml file are taken from the default values.
        The policies and stocks specified in the yaml file are directly loaded.

        :param filename: the path to a yaml file
        :type filename: str
        :return: a Configuration object with settings from the yaml file
        :rtype: Configuration
        """
        config_obj = Configuration()
        with open(filename, "r") as fileobj:
            _config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)

        config_obj._update_from_config(_config)

        for key, policy_spec in _config.get("policy", {}).get("files", {}).items():
            modelfile, templatefile = policy_spec
            config_obj.policy.load_policy(modelfile, templatefile, key)

        for key, stockfile in _config.get("stock", {}).get("files", {}).items():
            config_obj.stock.load_stock(stockfile, key)

        if "mongodb" in _config.get("stock", {}):
            query_obj = MongoDbInchiKeyQuery(**(_config["stock"]["mongodb"] or {}))
            config_obj.stock.load_stock(query_obj, "mongodb_stock")

        return config_obj

    def update(self, **settings):
        """ Update the configuration using dictionary of parameters
        """
        for setting, value in settings.items():
            setattr(self, setting, value)
            self._logger.info(f"Setting {setting.replace('_', ' ')} to {value}")

    def _update_from_config(self, config):
        self._properties.update(config.get("finder", {}).get("properties", {}))
        self._properties.update(config.get("policy", {}).get("properties", {}))
        self._properties.update(config.get("properties", {}))
        self.__dict__.update(self._properties)
