""" Module containing base classes used to score routes.
"""

from __future__ import annotations

import abc
import functools
from collections.abc import Sequence as SequenceAbc
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from rxnutils.routes.readers import read_aizynthfinder_dict

from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.utils.exceptions import ScorerException

if TYPE_CHECKING:
    from rxnutils.routes.base import SynthesisRoute

    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import (
        Optional,
        Sequence,
        StrDict,
        Tuple,
        TypeVar,
        Union,
    )

    _Scoreable = TypeVar("_Scoreable", MctsNode, ReactionTree)
    _Scoreables = Sequence[_Scoreable]
    _ScorerItemType = Union[_Scoreables, _Scoreable]


@dataclass
class SquashScaler:
    """
    Squash function loosely adapted from a sigmoid function with parameters
    to modify and offset the shape

    :param slope: the slope of the midpoint
    :param xoffset: the offset of the midpoint along the x-axis
    :param yoffset: the offset of the curve along the y-axis
    """

    slope: float
    xoffset: float
    yoffset: float

    def __call__(self, val: float) -> float:
        return 1 / (1 + np.exp(self.slope * -(val - self.xoffset))) - self.yoffset


@dataclass
class MinMaxScaler:
    """
    Scaling function that normalises the value between 0 - 1,
    the reverse variable controls the direction of scaling,
    reverse should set to be true for rewards that need to be minimised
    the scale_factor could be used to adjust the scores when they are too small or too big

    :param val: the value that is being scaled
    :param min_val: minimum val param val could take
    :param max_val: maximum val param val could take
    :param scale_factor: scaling factor applied to the minmax scaled output
    """

    min_val: float
    max_val: float
    reverse: bool
    scale_factor: float = 1

    def __call__(self, val: float) -> float:
        val = np.clip(val, self.min_val, self.max_val)
        if self.reverse:
            numerator = self.max_val - val
        else:
            numerator = val - self.min_val
        return (numerator / (self.max_val - self.min_val)) * self.scale_factor


@dataclass
class PowerScaler:
    """
    Scaling function that returns the base_coefficient to the power of the value.

    :param val: the value that is being scaled
    :param base_coefficient: coefficient for scalingm 'base_coefficient < 1' => reversing scale
    """

    base_coefficient: float = 0.98

    def __call__(self, val: float) -> float:
        return self.base_coefficient**val


_SCALERS = {"squash": SquashScaler, "min_max": MinMaxScaler, "power": PowerScaler}


class Scorer(abc.ABC):
    """
    Abstract base class for classes that do scoring on MCTS-like nodes or reaction trees.

    The actual scoring is done be calling an instance of
    a scorer class with a ``Node`` or ``ReactionTree`` object as only argument.

    .. code-block::

        scorer = MyScorer()
        score = scorer(node1)

    You can also give a list of such objects to the scorer

    .. code-block::

        scorer = MyScorer()
        scores = scorer([node1, node2])

    :param config: the configuration the tree search
    :param scaler_params: the parameter settings of the scaler
    """

    scorer_name = "base"

    def __init__(
        self,
        config: Optional[Configuration] = None,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        self._config = config
        self._reverse_order: bool = True
        self._scaler = None
        self._scaler_name = ""
        if scaler_params:
            self._scaler_name = scaler_params["name"]
            del scaler_params["name"]
            if scaler_params:
                self._scaler = _SCALERS[self._scaler_name](**scaler_params)
            else:
                # for parameterless function
                self._scaler = _SCALERS[self._scaler_name]()

    def __call__(self, item: _ScorerItemType) -> Union[float, Sequence[float]]:
        if isinstance(item, SequenceAbc):
            return self._score_many(item)
        if isinstance(item, (MctsNode, ReactionTree)):
            return self._score_just_one(item)  # type: ignore
        raise ScorerException(
            f"Unable to score item from class {item.__class__.__name__}"
        )

    def __repr__(self) -> str:
        repr_name = self.scorer_name
        if self._scaler_name:
            repr_name += f"-{self._scaler_name}"
        return repr_name

    def sort(
        self, items: _Scoreables
    ) -> Tuple[_Scoreables, Sequence[float], Sequence[int]]:
        """
        Sort nodes or reaction trees in descending order based on the score

        :param items: the items to sort
        :return: the sorted items and their scores
        """
        scores = self._score_many(items)
        assert isinstance(scores, SequenceAbc)
        sortidx = sorted(
            range(len(scores)), key=scores.__getitem__, reverse=self._reverse_order
        )
        scores = [scores[idx] for idx in sortidx]
        sorted_items = [items[idx] for idx in sortidx]
        return sorted_items, scores, sortidx

    def _score_just_one(self, item: _Scoreable) -> float:
        if isinstance(item, MctsNode):
            node_score = self._score_node(item)
            if self._scaler:
                node_score = self._scaler(node_score)
            return node_score
        if isinstance(item, ReactionTree):
            tree_score = self._score_reaction_tree(item)
            if self._scaler:
                tree_score = self._scaler(tree_score)
            return tree_score
        raise ScorerException(
            f"Unable to score item from class {item.__class__.__name__}"
        )

    def _score_many(self, items: _Scoreables) -> Sequence[float]:
        return [self._score_just_one(item) for item in items]

    @abc.abstractmethod
    def _score_node(self, node: MctsNode) -> float:
        pass

    @abc.abstractmethod
    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        pass


@functools.lru_cache
def make_rxnutils_route(tree: ReactionTree) -> SynthesisRoute:
    return read_aizynthfinder_dict(tree.to_dict())
