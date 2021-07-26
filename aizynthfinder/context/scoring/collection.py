""" Module containing classes used to score the reaction routes.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.context.scoring.scorers import __name__ as scorers_module
from aizynthfinder.context.scoring.scorers import (
    Scorer,
    StateScorer,
    NumberOfReactionsScorer,
    NumberOfPrecursorsScorer,
    NumberOfPrecursorsInStockScorer,
    AverageTemplateOccurrenceScorer,
)
from aizynthfinder.utils.exceptions import ScorerException
from aizynthfinder.utils.loading import load_dynamic_class

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Union,
        List,
        Any,
        Sequence,
        TypeVar,
    )
    from aizynthfinder.context.config import Configuration

    _Scoreable = TypeVar("_Scoreable", MctsNode, ReactionTree)
    _Scoreables = Sequence[_Scoreable]
    _ScorerItemType = Union[_Scoreables, _Scoreable]


_SIMPLE_SCORERS = [
    StateScorer,
    NumberOfReactionsScorer,
    NumberOfPrecursorsScorer,
    NumberOfPrecursorsInStockScorer,
    AverageTemplateOccurrenceScorer,
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
        for cls in _SIMPLE_SCORERS:
            self.load(cls(config))

    def load(self, scorer: Scorer) -> None:  # type: ignore
        """
        Add a pre-initialized scorer object to the collection

        :param scorer: the item to add
        """
        if not isinstance(scorer, Scorer):
            raise ScorerException(
                "Only objects of classes inherited from Scorer can be added"
            )
        self._items[repr(scorer)] = scorer

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

    def names(self) -> List[str]:
        """Return a list of the names of all the loaded scorers"""
        return self.items

    def objects(self) -> List[Scorer]:
        """Return a list of all the loaded scorer objects"""
        return list(self._items.values())
