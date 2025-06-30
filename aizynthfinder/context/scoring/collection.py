""" Module containing classes used to score the reaction routes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.context.scoring.scorers import StateScorer
from aizynthfinder.context.scoring.scorers import __name__ as scorers_module
from aizynthfinder.context.scoring.scorers_mols import (
    NumberOfPrecursorsInStockScorer,
    NumberOfPrecursorsScorer,
)
from aizynthfinder.context.scoring.scorers_reactions import (
    NumberOfReactionsScorer,
    Scorer,
)
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.utils.exceptions import ScorerException
from aizynthfinder.utils.loading import load_dynamic_class

if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Any, List, Sequence, TypeVar, Union

    _Scoreable = TypeVar("_Scoreable", MctsNode, ReactionTree)
    _Scoreables = Sequence[_Scoreable]
    _ScorerItemType = Union[_Scoreables, _Scoreable]


_SIMPLE_SCORERS = [
    StateScorer,
    NumberOfReactionsScorer,
    NumberOfPrecursorsScorer,
    NumberOfPrecursorsInStockScorer,
]


class ScorerCollection(ContextCollection):
    """
    Store scorer objects for the aizynthfinder interface.

    The scorers can be obtained by name with simple indexing

    .. code-block::

        scorers = ScorerCollection()
        scorer = scorers['state score']

    Scorers defined in this module and that does not require any
    other argument to initialize than the ``config`` are auto-loaded.

    :param config: the configuration of the tree search
    """

    _collection_name = "scorer"

    def __init__(self, config: Configuration) -> None:
        super().__init__()
        self._config = config
        self.create_default_scorers()

    def __repr__(self) -> str:
        if self.selection:
            return f"{self._collection_name} ({', '.join(self.selection)})"

        return f"{self._collection_name} ({', '.join(self.items)})"

    def create_default_scorers(self) -> None:
        """
        Setup the scores that only need the config as their input.
        """
        for cls in _SIMPLE_SCORERS:
            self.load(cls(self._config), silent=True)

    def load(self, scorer: Scorer, silent: bool = False) -> None:  # type: ignore
        """
        Add a pre-initialized scorer object to the collection

        :param scorer: the item to add
        :param silent: if True will not write out the name of the loaded scorer
        """
        if not isinstance(scorer, Scorer):
            raise ScorerException(
                "Only objects of classes inherited from Scorer can be added"
            )
        self._items[repr(scorer)] = scorer
        if not silent:
            self._logger.info(f"Loaded scorer: {repr(scorer)}")

    def load_from_config(self, **scorers_config: Any) -> None:
        """
        Load one or several scorers from a configuration dictionary

        The keys are the name of scorer class. If a scorer is not
        defined in the ``aizynthfinder.context.scoring`` module, the module
        name can be appended, e.g. ``mypackage.scoring.AwesomeScore``.

        The values of the configuration is passed directly to the scorer
        class along with the ``config`` parameter.

        :raises ScorerException: if module or class could not be found
        """
        for name_spec, scorer_config in scorers_config.items():
            config_str = (
                f" from configuration '{scorer_config}'" if scorer_config else ""
            )
            cls = load_dynamic_class(name_spec, scorers_module, ScorerException)
            obj = cls(self._config, **(scorer_config or {}))
            self._logger.info(f"Loaded scorer: '{repr(obj)}'{config_str}")
            self._items[repr(obj)] = obj

    def make_subset(self, subset_names: List[str]) -> "ScorerCollection":
        """
        Make a new scorer collection by taking a subset of this
        collection. The scorer instances will be shared between
        the collections

        :param subset_names: the scorers to copy over
        :returns: the newly formed collection
        """
        new_collection = ScorerCollection(self._config)
        for name in new_collection.names():
            del new_collection[name]
        for name in subset_names:
            try:
                scorer = self._items[name]
            except KeyError:
                raise ScorerException(f"Unable to find '{name}' in parent collection")
            new_collection.load(scorer, silent=True)
            new_collection.select(name, append=True)
        return new_collection

    def names(self) -> List[str]:
        """Return a list of the names of all the loaded scorers"""
        return self.items

    def objects(self) -> List[Scorer]:
        """Return a list of all the loaded scorer objects"""
        return list(self._items.values())

    def score_vector(self, item: _Scoreable) -> Sequence[float]:
        """
        For the given item, score it with all selected scorers
        and return a vector

        :param item: the item to be scored
        :returns: the vector with the scores
        """
        if not self.selection:
            return []
        return [self[scorer](item) for scorer in self.selection]

    def weighted_score(self, item: _Scoreable, weights: Sequence[float]) -> float:
        """
        For the given item, score it with all selected scorers
        and return a weighted sum of all the scores.

        If weights is not the same length as the number of scorers
        an exception is raised.

        If no scorers are selected this will raise an exception

        :param item: the item to be scored
        :param weights: the weights of the scorers
        :returns: the weighted sum
        """
        if not self.selection:
            raise ScorerException(
                "No scorers are selected so cannot compute weighted sum"
            )

        if len(weights) != len(self.selection):
            raise ScorerException(
                "The number of weights given does not agree with the number of scorers"
            )
        return sum(
            weight * score for weight, score in zip(weights, self.score_vector(item))
        )
