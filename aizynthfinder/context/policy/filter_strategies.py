""" Module containing classes that implements different filter policy strategies
"""
from __future__ import annotations
import abc
from typing import TYPE_CHECKING

import numpy as np

from aizynthfinder.utils.models import load_model
from aizynthfinder.utils.logging import logger
from aizynthfinder.context.policy.utils import _make_fingerprint
from aizynthfinder.utils.exceptions import PolicyException, RejectionException

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Any, List, Tuple
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem.reaction import RetroReaction


class FilterStrategy(abc.ABC):
    """
    A base class for all filter strategies.

    The filter can be applied by either calling the `apply` method
    of by calling the instantiated class with a reaction.

    .. code-block::

        filter = MyFilterStrategy("dummy", config)
        filter.apply(reaction)
        filter(reaction)

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    _required_kwargs: List[str] = []

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise PolicyException(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self._config = config
        self._logger = logger()
        self.key = key

    def __call__(self, reaction: RetroReaction) -> None:
        self.apply(reaction)

    @abc.abstractmethod
    def apply(self, reaction: RetroReaction) -> None:
        """
        Apply the filter on the reaction. If the reaction
        should be rejected a `RejectionException` is raised

        :param reaction: the reaction to filter
        :raises: if the reaction should be rejected.
        """


class QuickKerasFilter(FilterStrategy):
    """
    Filter quick-filter trained on artificial negative data

    :param key: the key or label
    :param config: the configuration of the tree search
    :param source: the source of the policy model
    """

    _required_kwargs: List[str] = ["source"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        super().__init__(key, config, **kwargs)
        source = kwargs["source"]
        self._logger.info(f"Loading filter policy model from {source} to {key}")
        self.model = load_model(source, key, self._config.use_remote_models)

    def apply(self, reaction: RetroReaction) -> None:
        feasible, prob = self.feasibility(reaction)
        if not feasible:
            raise RejectionException(f"{reaction} was filtered out with prob {prob}")

    def feasibility(self, reaction: RetroReaction) -> Tuple[bool, float]:
        """
        Computes if a given reaction is feasible by given
        the reaction fingerprint to a network model

        :param reaction: the reaction to query
        :return: if the reaction is feasible
        """
        if not reaction.reactants:
            return False, 0.0

        prob = self._predict(reaction)
        feasible = prob >= self._config.filter_cutoff
        return feasible, prob

    def _predict(self, reaction: RetroReaction) -> float:
        prod_fp, rxn_fp = self._reaction_to_fingerprint(reaction, self.model)
        return self.model.predict([prod_fp, rxn_fp])[0][0]

    @staticmethod
    def _reaction_to_fingerprint(
        reaction: RetroReaction, model: Any
    ) -> Tuple[np.ndarray, np.ndarray]:
        rxn_fp = _make_fingerprint(reaction, model)
        prod_fp = _make_fingerprint(reaction.mol, model)
        return prod_fp, rxn_fp
