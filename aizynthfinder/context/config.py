""" Module containing a class for encapsulating the settings of the tree search
"""
import os

import yaml

from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.paths import data_path
from aizynthfinder.context.policy import ExpansionPolicy, FilterPolicy
from aizynthfinder.context.stock import Stock
from aizynthfinder.context.scoring import ScorerCollection


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
        self.expansion_policy = ExpansionPolicy(self)
        self.filter_policy = FilterPolicy(self)
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

        config_obj.expansion_policy.load_from_config(
            **source.get("policy", {}).get("files", {})
        )
        config_obj.filter_policy.load_from_config(
            **source.get("filter", {}).get("files", {})
        )
        config_obj.stock.load_from_config(**source.get("stock", {}))

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
        """
        Update the configuration using dictionary of parameters

        If a value is None that setting is ignored.
        """
        for setting, value in settings.items():
            if value is None:
                continue
            setattr(self, setting, value)
            self._logger.info(f"Setting {setting.replace('_', ' ')} to {value}")

    def _update_from_config(self, config):
        self._properties.update(config.get("finder", {}).get("properties", {}))
        self._properties.update(config.get("policy", {}).get("properties", {}))
        self._properties.update(config.get("filter", {}).get("properties", {}))
        self._properties.update(config.get("properties", {}))
        self.__dict__.update(self._properties)
