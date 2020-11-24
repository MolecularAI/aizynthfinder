""" Module containing classes used to score the reaction routes.
"""
import abc
import importlib

from aizynthfinder.mcts.node import Node
from aizynthfinder.analysis import ReactionTree
from aizynthfinder.mcts.state import State
from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.collection import ContextCollection


class ScorerException(Exception):
    """ Exception raised by classes in this module
    """


class Scorer(abc.ABC):
    """
    Abstract base class for classes that do scoring on nodes or reaction trees.

    The actual scoring is done be calling an instance of
    a scorer class with a ``Node`` or ``ReactionTree`` object as only argument.

    .. code-block::

        scorer = MyScore()
        score = scorer(node1)

    :param config: the configuration the tree search
    :type abc: Configuration, optional
    """

    def __init__(self, config=None):
        self._config = config
        self._reverse_order = True

    def __call__(self, item):
        if isinstance(item, Node):
            return self._score_node(item)
        if isinstance(item, ReactionTree):
            return self._score_reaction_tree(item)
        raise ScorerException(
            f"Unable to score item from class {item.__class__.__name__}"
        )

    def sort(self, items, return_sort_indices=False):
        """
        Sort nodes or reaction trees in descending order based on the score

        :param items: the items to sort
        :type items: list of Nodes or list of ReactionTree
        :param return_sort_indices: if True, also returns the indices of the original list that sorts any list
        :type return_sort_indices: bool, optional
        :return: the sorted items and their scores
        :rtype: tuple of list of Nodes or list of ReactionTree, list of scores
        """
        scores = [self(item) for item in items]
        sortidx = sorted(
            range(len(scores)), key=scores.__getitem__, reverse=self._reverse_order
        )
        scores = [scores[idx] for idx in sortidx]
        sorted_items = [items[idx] for idx in sortidx]
        if return_sort_indices:
            return sorted_items, scores, sortidx
        return sorted_items, scores

    @abc.abstractmethod
    def _score_node(self, node):
        pass

    @abc.abstractmethod
    def _score_reaction_tree(self, tree):
        pass


class StateScorer(Scorer):
    """ Class for scoring nodes based on the state score
    """

    def __repr__(self):
        return "state score"

    def _score_node(self, node):
        return node.state.score

    def _score_reaction_tree(self, tree):
        mols = [
            TreeMolecule(
                parent=None, transform=tree.depth(leaf) // 2, smiles=leaf.smiles
            )
            for leaf in tree.leafs()
        ]
        state = State(mols, self._config)
        return state.score


class NumberOfReactionsScorer(Scorer):
    """ Class for scoring nodes based on the number of reaction it took to get to a node
    """

    def __init__(self, config=None):
        super().__init__(config)
        self._reverse_order = False

    def __repr__(self):
        return "number of reactions"

    def _score_node(self, node):
        reactions, _ = node.tree.route_to_node(node)
        return len(reactions)

    def _score_reaction_tree(self, tree):
        return len(list(tree.reactions()))


class NumberOfPrecursorsScorer(Scorer):
    """ Class for scoring nodes based on the number of pre-cursors in a node or route
    """

    def __init__(self, config=None):
        super().__init__(config)
        self._reverse_order = False

    def __repr__(self):
        return "number of pre-cursors"

    def _score_node(self, node):
        return len(node.state.mols)

    def _score_reaction_tree(self, tree):
        return len(list(tree.leafs()))


class NumberOfPrecursorsInStockScorer(Scorer):
    """ Class for scoring nodes based on the number of pre-cursors in stock a node or route
    """

    def __repr__(self):
        return "number of pre-cursors in stock"

    def _score_node(self, node):
        return len([mol for mol in node.state.mols if mol in self._config.stock])

    def _score_reaction_tree(self, tree):
        return len([mol for mol in tree.leafs() if mol in self._config.stock])


class AverageTemplateOccurenceScorer(Scorer):
    """ Class for scoring the nodes based on the average occurence of the templates used to get to a node
    """

    def __repr__(self):
        return "average template occurence"

    def _score_node(self, node):
        reactions, _ = node.tree.route_to_node(node)
        if not reactions:
            return 0.0
        occurences = [
            reaction.metadata.get("library_occurence", 0) for reaction in reactions
        ]
        return sum(occurences) / len(reactions)

    def _score_reaction_tree(self, tree):
        reactions = list(tree.reactions())
        if not reactions:
            return 0.0
        occurences = [
            reaction.metadata.get("library_occurence", 0) for reaction in reactions
        ]
        return sum(occurences) / len(reactions)


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
    :type config: Configuration
    """

    _collection_name = "scorer"

    def __init__(self, config):
        super().__init__()
        self._config = config
        for cls in _SIMPLE_SCORERS:
            self.load(cls(config))

    def load(self, scorer):
        """
        Add a pre-initialized scorer object to the collection

        :param scorer: the item to add
        :type scorer: Scorer
        """
        if not isinstance(scorer, Scorer):
            raise ScorerException(
                "Only objects of classes inherited from Scorer can be added"
            )
        self._items[repr(scorer)] = scorer

    def load_from_config(self, **scorers_config):
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

            obj = getattr(loaded_module, name)(self._config, **scorer_config)
            self._items[repr(obj)] = obj

    def names(self):
        """ Return a list of the names of all the loaded scorers
        """
        return self.items

    def objects(self):
        """ Return a list of all the loaded scorer objects
        """
        return list(self._items.values())
