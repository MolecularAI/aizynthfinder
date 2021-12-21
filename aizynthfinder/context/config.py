""" Module containing a class for encapsulating the settings of the tree search
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING

import yaml

from aizynthfinder.utils.logging import logger
from aizynthfinder.context.policy import ExpansionPolicy, FilterPolicy
from aizynthfinder.context.stock import Stock
from aizynthfinder.context.scoring import ScorerCollection
from aizynthfinder.context.cost import MoleculeCost


if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import StrDict, Any, Dict, Union, Optional


@dataclass
class _PostprocessingConfiguration:
    min_routes: int = 5
    max_routes: int = 25
    all_routes: bool = False
    route_distance_model: Optional[str] = None


@dataclass
class Configuration:  # pylint: disable=R0902
    """
    Encapsulating the settings of the tree search, including the policy,
    the stock, the loaded scorers and various parameters.

    All the parameters can be retrieved as attributes of the Configuration
    object, e.g.

    .. code-block::

        config.max_transforms  # The maximum of transform
        config.iteration_limit  # The maximum number of iterations
    """

    C: float = 1.4  # pylint: disable=invalid-name
    cutoff_cumulative: float = 0.995
    cutoff_number: int = 50
    additive_expansion: bool = False
    use_rdchiral: bool = True
    max_transforms: int = 6
    default_prior: float = 0.5
    use_prior: bool = True
    iteration_limit: int = 100
    return_first: bool = False
    time_limit: int = 120
    filter_cutoff: float = 0.05
    exclude_target_from_stock: bool = True
    template_column: str = "retro_template"
    prune_cycles_in_search: bool = True
    use_remote_models: bool = False
    search_algorithm: str = "mcts"
    post_processing: _PostprocessingConfiguration = _PostprocessingConfiguration()
    stock: Stock = None  # type: ignore
    expansion_policy: ExpansionPolicy = None  # type: ignore
    filter_policy: FilterPolicy = None  # type: ignore
    scorers: ScorerCollection = None  # type: ignore
    molecule_cost: MoleculeCost = None  # type: ignore

    def __post_init__(self) -> None:
        self._properties: StrDict = {}
        self.stock = Stock()
        self.expansion_policy = ExpansionPolicy(self)
        self.filter_policy = FilterPolicy(self)
        self.scorers = ScorerCollection(self)
        self.molecule_cost = MoleculeCost()
        self._logger = logger()

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Configuration):
            return False
        return self.properties == other.properties

    @classmethod
    def from_dict(cls, source: StrDict) -> "Configuration":
        """
        Loads a configuration from a dictionary structure.
        The parameters not set in the dictionary are taken from the default values.
        The policies and stocks specified are directly loaded.

        :param source: the dictionary source
        :return: a Configuration object with settings from the source
        """
        # pylint: disable=protected-access
        config_obj = Configuration()
        src_copy = dict(source)  # Make a copy so we can pop items
        config_obj._update_from_config(src_copy)

        config_obj.expansion_policy.load_from_config(**src_copy.get("policy", {}))
        config_obj.filter_policy.load_from_config(**src_copy.get("filter", {}))
        config_obj.stock.load_from_config(**src_copy.get("stock", {}))
        config_obj.scorers.load_from_config(**src_copy.get("scorer", {}))
        config_obj.molecule_cost.load_from_config(**src_copy.get("molecule_cost", {}))

        return config_obj

    @classmethod
    def from_file(cls, filename: str) -> "Configuration":
        """
        Loads a configuration from a yaml file.
        The parameters not set in the yaml file are taken from the default values.
        The policies and stocks specified in the yaml file are directly loaded.

        :param filename: the path to a yaml file
        :return: a Configuration object with settings from the yaml file
        """
        with open(filename, "r") as fileobj:
            _config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
        return Configuration.from_dict(_config)

    @property
    def properties(self) -> Dict[str, Union[int, float, str, bool]]:
        """Return the basic properties of the config as a dictionary"""
        dict_ = {}
        for item in dir(self):
            if item == "properties" or item.startswith("_"):
                continue
            attr = getattr(self, item)
            if isinstance(attr, (int, float, str, bool)):
                dict_[item] = attr
        return dict_

    @properties.setter
    def properties(self, dict_: Dict[str, Union[int, float, str, bool]]) -> None:
        """
        Update the configuration using dictionary of parameters

        If a value is None that setting is ignored.

        :param dict_: the dictionary of properties
        """
        for setting, value in dict_.items():
            if value is None:
                continue
            if not hasattr(self, setting):
                raise AttributeError(f"Could not find attribute to set: {setting}")
            if not isinstance(value, (int, float, str, bool)):
                raise ValueError(
                    f"Trying to set property with invalid value {setting}={value}"
                )
            setattr(self, setting, value)
            self._logger.info(f"Setting {setting.replace('_', ' ')} to {value}")

    def _update_from_config(self, config: StrDict) -> None:
        #  The first 3 places for properties are kept for historical reasons, but they are not recommended usage
        dict_ = config.get("finder", {}).pop("properties", {})
        dict_.update(config.get("policy", {}).pop("properties", {}))
        dict_.update(config.get("filter", {}).pop("properties", {}))
        dict_.update(config.pop("properties", {}))
        self.post_processing = _PostprocessingConfiguration(
            **dict_.pop("post_processing", {})
        )
        self.properties = dict_
