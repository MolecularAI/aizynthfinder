""" Module containing a class for encapsulating the settings of the tree search
"""
import os
import importlib

import yaml

from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.paths import data_path
from aizynthfinder.mcts.policy import Policy
from aizynthfinder.mcts.stock import Stock, MongoDbInchiKeyQuery
from aizynthfinder.scoring import ScorerCollection


class Configuration:
    """
    Encapsulating the settings of the tree search, including the policy,
    the stock, the loaded scorers and various parameters.

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
        self.scorers = ScorerCollection(self)
        self._logger = logger()

    def __eq__(self, other):
        return self._properties == other._properties

    @classmethod
    def from_dict(cls, source):
        """
        Loads a configuration from a dictionary structure.
        The parameters not set in the dictionary are taken from the default values.
        The policies and stocks specified are directly loaded.

        :param source: the dictionary source
        :type source: dict
        :return: a Configuration object with settings from the source
        :rtype: Configuration
        """
        config_obj = Configuration()
        config_obj._update_from_config(source)

        for key, policy_spec in source.get("policy", {}).get("files", {}).items():
            modelfile, templatefile = policy_spec
            config_obj.policy.load_policy(modelfile, templatefile, key)

        config_obj._load_stocks(source)

        return config_obj

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
        with open(filename, "r") as fileobj:
            _config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
        return Configuration.from_dict(_config)

    def update(self, **settings):
        """ Update the configuration using dictionary of parameters
        """
        for setting, value in settings.items():
            setattr(self, setting, value)
            self._logger.info(f"Setting {setting.replace('_', ' ')} to {value}")

    def _load_stocks(self, config):
        for key, stockfile in config.get("stock", {}).get("files", {}).items():
            self.stock.load_stock(stockfile, key)

        if "mongodb" in config.get("stock", {}):
            query_obj = MongoDbInchiKeyQuery(**(config["stock"]["mongodb"] or {}))
            self.stock.load_stock(query_obj, "mongodb_stock")

        # Load stocks specifying a module and class, e.g. package.module.MyQueryClass
        for name, stock_config in config.get("stock", {}).items():
            if name in ["files", "mongodb"] or name.find(".") == -1:
                continue

            module_name, class_name = name.rsplit(".", maxsplit=1)
            try:
                module = importlib.import_module(module_name)
            except ImportError:
                self._logger.warning(
                    f"Unable to load module '{module_name}' containing stock query classes"
                )
                continue

            if hasattr(module, class_name):
                query_obj = getattr(module, class_name)(**(stock_config or {}))
                self.stock.load_stock(query_obj, class_name)
            else:
                self._logger.warning(
                    f"Unable to find query class '{class_name}' in '{module_name}''"
                )

    def _update_from_config(self, config):
        self._properties.update(config.get("finder", {}).get("properties", {}))
        self._properties.update(config.get("policy", {}).get("properties", {}))
        self._properties.update(config.get("properties", {}))
        self.__dict__.update(self._properties)
