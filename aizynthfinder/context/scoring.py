""" Module containing classes used to score the reaction routes.
"""
from __future__ import annotations
import abc
import importlib
from collections import defaultdict
from collections.abc import Sequence as SequenceAbc
from typing import TYPE_CHECKING

from aizynthfinder.mcts.node import Node
from aizynthfinder.analysis import ReactionTree
from aizynthfinder.mcts.state import State
from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.collection import ContextCollection
from aizynthfinder.context.stock import StockException

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Optional,
        Union,
        List,
        Tuple,
        Any,
        Sequence,
        Iterable,
        TypeVar,
    )
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem import Molecule

    _Scoreable = TypeVar("_Scoreable", Node, ReactionTree)
    _Scoreables = Sequence[_Scoreable]
    _ScorerItemType = Union[_Scoreables, _Scoreable]


class ScorerException(Exception):
    """Exception raised by classes in this module"""


class Scorer(abc.ABC):
    """
    Abstract base class for classes that do scoring on nodes or reaction trees.

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
    """

    def __init__(self, config: Configuration = None) -> None:
        self._config = config
        self._reverse_order: bool = True

    def __call__(self, item: _ScorerItemType) -> Union[float, Sequence[float]]:
        if isinstance(item, SequenceAbc):
            return self._score_many(item)
        if isinstance(item, Node) or isinstance(item, ReactionTree):
            return self._score_just_one(item)  # type: ignore
        raise ScorerException(
            f"Unable to score item from class {item.__class__.__name__}"
        )

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
        if isinstance(item, Node):
            return self._score_node(item)
        return self._score_reaction_tree(item)

    def _score_many(self, items: _Scoreables) -> Sequence[float]:
        if isinstance(items[0], Node):
            return self._score_nodes(items)
        elif isinstance(items[0], ReactionTree):
            return self._score_reaction_trees(items)
        raise ScorerException(
            f"Unable to score item from class {items[0].__class__.__name__}"
        )

    @abc.abstractmethod
    def _score_node(self, node: Node) -> float:
        pass

    def _score_nodes(self, nodes: _Scoreables) -> Sequence[float]:
        return [self._score_node(node) for node in nodes]

    @abc.abstractmethod
    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        pass

    def _score_reaction_trees(self, trees: _Scoreables) -> Sequence[float]:
        return [self._score_reaction_tree(tree) for tree in trees]


class StateScorer(Scorer):
    """Class for scoring nodes based on the state score"""

    def __init__(self, config: Configuration) -> None:
        super().__init__(config)
        self._config: Configuration = config

    def __repr__(self) -> str:
        return "state score"

    def _score_node(self, node: Node) -> float:
        return node.state.score

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        mols = [
            TreeMolecule(
                parent=None, transform=tree.depth(leaf) // 2, smiles=leaf.smiles
            )
            for leaf in tree.leafs()
        ]
        state = State(mols, self._config)
        return state.score


class NumberOfReactionsScorer(Scorer):
    """Class for scoring nodes based on the number of reaction it took to get to a node"""

    def __init__(self, config: Configuration = None) -> None:
        super().__init__(config)
        self._reverse_order = False

    def __repr__(self) -> str:
        return "number of reactions"

    def _score_node(self, node: Node) -> float:
        reactions, _ = node.tree.route_to_node(node)
        return len(reactions)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.reactions()))


class NumberOfPrecursorsScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in a node or route"""

    def __init__(self, config: Configuration = None) -> None:
        super().__init__(config)
        self._reverse_order = False

    def __repr__(self) -> str:
        return "number of pre-cursors"

    def _score_node(self, node: Node) -> float:
        return len(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.leafs()))


class NumberOfPrecursorsInStockScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in stock a node or route"""

    def __init__(self, config: Configuration) -> None:
        super().__init__(config)
        self._stock = config.stock

    def __repr__(self) -> str:
        return "number of pre-cursors in stock"

    def _score_node(self, node: Node) -> float:
        return len([mol for mol in node.state.mols if mol in self._stock])

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len([mol for mol in tree.leafs() if mol in self._stock])


class AverageTemplateOccurenceScorer(Scorer):
    """Class for scoring the nodes based on the average occurrence of the templates used to get to a node"""

    def __repr__(self) -> str:
        return "average template occurence"

    def _score_node(self, node: Node) -> float:
        reactions, _ = node.tree.route_to_node(node)
        if not reactions:
            return 0.0
        occurences = [
            reaction.metadata.get("library_occurence", 0) for reaction in reactions
        ]
        return sum(occurences) / len(reactions)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        reactions = list(tree.reactions())
        if not reactions:
            return 0.0
        occurences = [
            reaction.metadata.get("library_occurence", 0) for reaction in reactions
        ]
        return sum(occurences) / len(reactions)


class PriceSumScorer(Scorer):
    def __init__(
        self,
        config: Configuration,
        default_cost: float = 1.0,
        not_in_stock_multiplier: int = 10,
    ) -> None:
        super().__init__(config)
        self._config: Configuration = config

        self.default_cost = default_cost
        self.not_in_stock_multiplier = not_in_stock_multiplier
        self._reverse_order = False

    def __repr__(self) -> str:
        return "sum of prices"

    def _calculate_leaf_costs(
        self, leafs: Union[Sequence[Molecule], Iterable[Molecule]]
    ) -> dict:
        costs = {}
        for mol in leafs:
            if mol not in self._config.stock:
                continue
            try:
                cost = self._config.stock.price(mol)
            except StockException:
                costs[mol] = self.default_cost
            else:
                costs[mol] = cost

        max_cost = max(costs.values()) if costs else self.default_cost
        return defaultdict(lambda: max_cost * self.not_in_stock_multiplier, costs)

    def _score_node(self, node: Node) -> float:
        leaf_costs = self._calculate_leaf_costs(node.state.mols)
        return sum(leaf_costs[mol] for mol in node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        leaf_costs = self._calculate_leaf_costs(tree.leafs())
        return sum(leaf_costs[leaf] for leaf in tree.leafs())


class RouteCostScorer(PriceSumScorer):
    """
    Score based on the cost of molecules and reactions.
    From Badowski et al. Chem Sci. 2019, 10, 4640
    """

    def __init__(
        self,
        config: Configuration,
        reaction_cost: int = 1,
        average_yield: float = 0.8,
        default_cost: int = 1,
        not_in_stock_multiplier: int = 10,
    ) -> None:
        super().__init__(
            config,
            default_cost=default_cost,
            not_in_stock_multiplier=not_in_stock_multiplier,
        )
        self.reaction_cost = reaction_cost
        self.average_yield = average_yield
        self._reverse_order = False

    def __repr__(self) -> str:
        return "route cost"

    def _score_node(self, node: Node) -> float:
        leaf_costs = self._calculate_leaf_costs(node.state.mols)

        reactions, nodes = node.tree.route_to_node(node)
        if not reactions:
            return leaf_costs[node.state.mols[0]]

        scores = {id(mol): leaf_costs[mol] for mol in nodes[-1].state.mols}
        for node, reaction in zip(nodes[::-1][1:], reactions[::-1]):
            updated_scores = {
                id(mol): scores[id(mol)]
                for mol in node.state.mols
                if mol != reaction.mol
            }
            child_sum = sum(
                1 / self.average_yield * score
                for id_, score in scores.items()
                if id_ not in updated_scores
            )
            updated_scores[id(reaction.mol)] = self.reaction_cost + child_sum
            scores = updated_scores

        return list(scores.values())[0]

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        def _recursive_score(node):
            # This list should contains 0 or 1 elements
            reaction_nodes = list(tree.graph[node])
            if not reaction_nodes:
                return leaf_costs[node]

            child_sum = sum(
                1 / self.average_yield * _recursive_score(child)
                for child in tree.graph[reaction_nodes[0]]
            )
            return self.reaction_cost + child_sum

        leaf_costs = self._calculate_leaf_costs(tree.leafs())
        return _recursive_score(tree.root)


_SIMPLE_SCORERS = [
    StateScorer,
    NumberOfReactionsScorer,
    NumberOfPrecursorsScorer,
    NumberOfPrecursorsInStockScorer,
    AverageTemplateOccurenceScorer,
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
            if "." not in name_spec:
                name = name_spec
                module_name = self.__module__
            else:
                module_name, name = name_spec.rsplit(".", maxsplit=1)

            try:
                loaded_module = importlib.import_module(module_name)
            except ImportError:
                raise ScorerException(f"Unable to load module: {module_name}")

            if not hasattr(loaded_module, name):
                raise ScorerException(
                    f"Module ({module_name}) does not have a class called {name}"
                )

            config_str = (
                f" from configuration '{scorer_config}'" if scorer_config else ""
            )
            obj = getattr(loaded_module, name)(self._config, **(scorer_config or {}))
            self._logger.info(f"Loaded scorer: '{repr(obj)}'{config_str}")
            self._items[repr(obj)] = obj

    def names(self) -> List[str]:
        """Return a list of the names of all the loaded scorers"""
        return self.items

    def objects(self) -> List[Scorer]:
        """Return a list of all the loaded scorer objects"""
        return list(self._items.values())
