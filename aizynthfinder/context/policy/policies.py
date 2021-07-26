""" Module containing classes that interfaces neural network policies
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.utils.loading import load_dynamic_class
from aizynthfinder.utils.exceptions import PolicyException
from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.context.policy.expansion_strategies import (
    ExpansionStrategy,
    TemplateBasedExpansionStrategy,
)
from aizynthfinder.context.policy.filter_strategies import (
    FilterStrategy,
    QuickKerasFilter,
)
from aizynthfinder.context.policy.expansion_strategies import (
    __name__ as expansion_strategy_module,
)
from aizynthfinder.context.policy.filter_strategies import (
    __name__ as filter_strategy_module,
)

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Any, Sequence, List, Tuple
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction


class ExpansionPolicy(ContextCollection):
    """
    An abstraction of an expansion policy.

    This policy provides actions that can be applied to a molecule

    :param config: the configuration of the tree search
    """

    _collection_name = "expansion policy"

    def __init__(self, config: Configuration) -> None:
        super().__init__()
        self._config = config

    def __call__(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[RetroReaction], List[float]]:
        return self.get_actions(molecules)

    def get_actions(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies

        :param molecules: the molecules to consider
        :return: the actions and the priors of those actions
        :raises: PolicyException: if the policy isn't selected
        """
        if not self.selection:
            raise PolicyException("No expansion policy selected")

        all_possible_actions = []
        all_priors = []
        for name in self.selection:
            possible_actions, priors = self[name].get_actions(molecules)
            all_possible_actions.extend(possible_actions)
            all_priors.extend(priors)
            if not self._config.additive_expansion and all_possible_actions:
                break
        return all_possible_actions, all_priors

    def load(self, source: ExpansionStrategy) -> None:  # type: ignore
        """
        Add a pre-initialized expansion strategy object to the policy

        :param source: the item to add
        """
        if not isinstance(source, ExpansionStrategy):
            raise PolicyException(
                "Only objects of classes inherited from ExpansionStrategy can be added"
            )
        self._items[source.key] = source

    def load_from_config(self, **config: Any) -> None:
        """
        Load one or more expansion policy from a configuration

        The format should be
        files:
            key:
                - path_to_model
                - path_to_templates
        or
        template-based:
            key:
                - path_to_model
                - path_to_templates
        or
        custom_package.custom_model.CustomClass:
            key:
                param1: value1
                param2: value2

        :param config: the configuration
        """
        files_spec = config.get("files", config.get("template-based", {}))
        for key, policy_spec in files_spec.items():
            modelfile, templatefile = policy_spec
            strategy = TemplateBasedExpansionStrategy(
                key, self._config, source=modelfile, templatefile=templatefile
            )
            self.load(strategy)

        # Load policies specifying a module and class, e.g. package.module.MyStrategyClass
        for strategy_spec, strategy_config in config.items():
            if strategy_spec in ["files", "template-based"]:
                continue
            cls = load_dynamic_class(
                strategy_spec, expansion_strategy_module, PolicyException
            )
            for key, policy_spec in strategy_config.items():
                obj = cls(key, self._config, **(policy_spec or {}))
                self.load(obj)


class FilterPolicy(ContextCollection):
    """
    An abstraction of a filter policy.

    This policy provides a query on a reaction to determine whether it should be rejected

    :param config: the configuration of the tree search
    """

    _collection_name = "filter policy"

    def __init__(self, config: Configuration) -> None:
        super().__init__()
        self._config = config

    def __call__(self, reaction: RetroReaction) -> None:
        return self.apply(reaction)

    def apply(self, reaction: RetroReaction) -> None:
        """
        Apply the all the selected filters on the reaction. If the reaction
        should be rejected a `RejectionException` is raised

        :param reaction: the reaction to filter
        :raises: if the reaction should be rejected or if a policy is selected
        """
        if not self.selection:
            raise PolicyException("No filter policy selected")

        for name in self.selection:
            self[name](reaction)

    def load(self, source: FilterStrategy) -> None:  # type: ignore
        """
        Add a pre-initialized filter strategy object to the policy

        :param source: the item to add
        """
        if not isinstance(source, FilterStrategy):
            raise PolicyException(
                "Only objects of classes inherited from FilterStrategy can be added"
            )
        self._items[source.key] = source

    def load_from_config(self, **config: Any) -> None:
        """
        Load one or more filter policy from a configuration

        The format should be
        files:
            key: path_to_model
        or
        quick-filter:
            key: path_to_model
        or
        custom_package.custom_model.CustomClass:
            key:
                param1: value1
                param2: value2

        :param config: the configuration
        """
        files_spec = config.get("files", config.get("quick-filter", {}))
        for key, modelfile in files_spec.items():
            strategy = QuickKerasFilter(key, self._config, source=modelfile)
            self.load(strategy)

        # Load policies specifying a module and class, e.g. package.module.MyStrategyClass
        for strategy_spec, strategy_config in config.items():
            if strategy_spec in ["files", "quick-filter"]:
                continue
            cls = load_dynamic_class(
                strategy_spec, filter_strategy_module, PolicyException
            )
            for key, policy_spec in strategy_config.items():
                obj = cls(key, self._config, **(policy_spec or {}))
                self.load(obj)
