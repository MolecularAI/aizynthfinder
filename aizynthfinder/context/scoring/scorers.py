""" Module containing classes used to score the reaction routes.
"""

from __future__ import annotations

import abc
import json
from collections import defaultdict
from collections.abc import Sequence as SequenceAbc
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

try:
    from route_distances.route_distances import route_distances_calculator
except ImportError:
    SUPPORT_DISTANCES = False
else:
    SUPPORT_DISTANCES = True

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.stock import StockException
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.utils.bonds import BrokenBonds
from aizynthfinder.utils.exceptions import ScorerException
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.sc_score import SCScore

if TYPE_CHECKING:
    from aizynthfinder.chem import (
        FixedRetroReaction,
        Molecule,
        RetroReaction,
        UniqueMolecule,
    )
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import (
        Iterable,
        Optional,
        Sequence,
        StrDict,
        Tuple,
        TypeVar,
        Union,
    )

    _Molecules = Sequence[Molecule]
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


_SCALERS = {"squash": SquashScaler, "min_max": MinMaxScaler}


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
                # for paramterless function
                self._scaler = _SCALERS[self._scaler_name]

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


class StateScorer(Scorer):
    """Class for scoring nodes based on the state score"""

    scorer_name = "state score"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        # This is necessary because config should not be optional for this scorer
        self._config: Configuration = config
        self._transform_scorer = MaxTransformScorerer(
            config,
            scaler_params={"name": "squash", "slope": -1, "yoffset": 0, "xoffset": 4},
        )
        self._in_stock_scorer = FractionInStockScorer(config)

    def _score(self, item: _Scoreable) -> float:
        in_stock_fraction = self._in_stock_scorer(item)
        max_transform = self._transform_scorer(item)
        # A scorer can return a list of float if the item is a list of trees/nodes,
        # but that is not the case here. However this is needed because of mypy
        assert isinstance(in_stock_fraction, float) and isinstance(max_transform, float)
        return 0.95 * in_stock_fraction + 0.05 * max_transform

    def _score_node(self, node: MctsNode) -> float:
        return self._score(node)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._score(tree)


class MaxTransformScorerer(Scorer):
    """Class for scoring nodes based on the maximum transform"""

    scorer_name = "max transform"

    def _score_node(self, node: MctsNode) -> float:
        return node.state.max_transforms

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        mols = [
            TreeMolecule(
                parent=None, transform=tree.depth(leaf) // 2, smiles=leaf.smiles
            )
            for leaf in tree.leafs()
        ]
        return max(mol.transform for mol in mols)


class FractionInStockScorer(Scorer):
    """Class for scoring nodes based on the fraction in stock"""

    scorer_name = "fraction in stock"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        # This is necessary because config should not be optional for this scorer
        self._config: Configuration = config

    def _score_node(self, node: MctsNode) -> float:
        num_in_stock = np.sum(node.state.in_stock_list)
        num_molecules = len(node.state.mols)
        return float(num_in_stock) / float(num_molecules)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        leaves = list(tree.leafs())
        num_in_stock = sum(mol in self._config.stock for mol in leaves)
        num_molecules = len(leaves)
        return float(num_in_stock) / float(num_molecules)


class NumberOfReactionsScorer(Scorer):
    """Class for scoring nodes based on the number of reaction it took to get to a node"""

    scorer_name = "number of reactions"

    def __init__(
        self,
        config: Optional[Configuration] = None,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._reverse_order = False

    def _score_node(self, node: MctsNode) -> float:
        reactions = node.actions_to()
        return len(reactions)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.reactions()))


class NumberOfPrecursorsScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in a node or route"""

    scorer_name = "number of pre-cursors"

    def __init__(
        self,
        config: Optional[Configuration] = None,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._reverse_order = False

    def _score_node(self, node: MctsNode) -> float:
        return len(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.leafs()))


class NumberOfPrecursorsInStockScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in stock a
    node or route"""

    scorer_name = "number of pre-cursors in stock"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        self._stock = config.stock

    def _score_node(self, node: MctsNode) -> float:
        return len([mol for mol in node.state.mols if mol in self._stock])

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len([mol for mol in tree.leafs() if mol in self._stock])


class AverageTemplateOccurrenceScorer(Scorer):
    """Class for scoring the nodes based on the average occurrence of the
    templates used to get to a node"""

    scorer_name = "average template occurrence"

    def _calc_average(
        self, reactions: Sequence[Union[FixedRetroReaction, RetroReaction]]
    ) -> float:
        if not reactions:
            return 0.0
        occurrences = [self._get_occurrence(reaction) for reaction in reactions]
        return sum(occurrences) / len(reactions)

    def _score_node(self, node: MctsNode) -> float:
        return self._calc_average(node.actions_to())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calc_average(list(tree.reactions()))

    @staticmethod
    def _get_occurrence(reaction: Union[FixedRetroReaction, RetroReaction]) -> int:
        return reaction.metadata.get(
            "library_occurrence", reaction.metadata.get("library_occurence", 0)
        )


class PriceSumScorer(Scorer):
    """Scorer that sums the prices of all pre-cursors"""

    scorer_name = "sum of prices"

    def __init__(
        self,
        config: Configuration,
        scaler_params: Optional[StrDict] = None,
        default_cost: float = 1.0,
        not_in_stock_multiplier: int = 10,
    ) -> None:
        super().__init__(config, scaler_params)
        self._config: Configuration = config
        self.default_cost = default_cost
        self.not_in_stock_multiplier = not_in_stock_multiplier
        self._reverse_order = False

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

    def _score_node(self, node: MctsNode) -> float:
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

    scorer_name = "route cost"

    def __init__(
        self,
        config: Configuration,
        scaler_params: Optional[StrDict] = None,
        reaction_cost: int = 1,
        average_yield: float = 0.8,
        default_cost: int = 1,
        not_in_stock_multiplier: int = 10,
    ) -> None:
        super().__init__(
            config,
            scaler_params,
            default_cost=default_cost,
            not_in_stock_multiplier=not_in_stock_multiplier,
        )
        self.reaction_cost = reaction_cost
        self.average_yield = average_yield
        self._reverse_order = False

    def _score_node(self, node: MctsNode) -> float:
        leaf_costs = self._calculate_leaf_costs(node.state.mols)

        reactions, nodes = node.path_to()
        if not reactions:
            return leaf_costs[node.state.mols[0]]

        scores = {id(mol): leaf_costs[mol] for mol in nodes[-1].state.mols}
        for pnode, reaction in zip(nodes[::-1][1:], reactions[::-1]):
            updated_scores = {
                id(mol): scores[id(mol)]
                for mol in pnode.state.mols
                if mol is not reaction.mol
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


class ReactionClassMembershipScorer(Scorer):
    """
    Scorer that checks if the reaction classes are in a specified set

    The score is calculated as product over each reaction. For each reaction
    the reaction classification is checked if it is in a given list of classes.
    """

    def __init__(
        self,
        config: Configuration,
        reaction_class_set: Sequence[str],
        in_set_score: float = 1.0,
        not_in_set_score: float = 0.1,
    ) -> None:
        super().__init__(config)
        self.reaction_class_set = set(item.split(" ")[0] for item in reaction_class_set)
        self.in_set_score = in_set_score
        self.not_in_set_score = not_in_set_score

    def __repr__(self) -> str:
        return "reaction class membership"

    def _calc_product(
        self, reactions: Sequence[Union[FixedRetroReaction, RetroReaction]]
    ) -> float:
        if not reactions:
            return 1.0
        prod = 1.0
        for reaction in reactions:
            membership = self._get_membership(reaction)
            prod *= self.in_set_score if membership else self.not_in_set_score
        return prod

    def _get_membership(
        self, reaction: Union[FixedRetroReaction, RetroReaction]
    ) -> bool:
        classification = reaction.metadata.get("classification")
        if not classification:
            return 0.0 in self.reaction_class_set
        return classification.split(" ")[0] in self.reaction_class_set

    def _score_node(self, node: MctsNode) -> float:
        return self._calc_product(node.actions_to())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calc_product(list(tree.reactions()))


class StockAvailabilityScorer(Scorer):
    """
    Scorer that computes score based on the stock availability of the starting material

    The score is calculated as a product of a stock score per starting material. The stock
    score for each molecule is based on the source of the stock, or a default value if the
    molecule is not in stock. The `other_source_score` parameter can be used to
    distinguish between "not in stock" and "not in the specificed sources" cases.
    """

    def __init__(
        self,
        config: Configuration,
        source_score: StrDict,
        default_score: float = 0.1,
        other_source_score: Optional[float] = None,
    ) -> None:
        super().__init__(config)
        assert self._config is not None
        self.source_score = source_score
        self.default_score = default_score
        self.other_source_score = other_source_score

    def __repr__(self) -> str:
        return "stock availability"

    def _calculate_leaf_costs(
        self, leafs: Union[Sequence[Molecule], Iterable[Molecule]]
    ) -> float:
        assert self._config is not None
        prod = 1.0
        for mol in leafs:
            availability = self._config.stock.availability_list(mol)
            scores = [
                self.source_score[source]
                for source in availability
                if source in self.source_score
            ]
            if scores:
                prod *= max(scores)
            elif self.other_source_score and mol in self._config.stock:
                prod *= self.other_source_score
            else:
                prod *= self.default_score
        return prod

    def _score_node(self, node: MctsNode) -> float:
        return self._calculate_leaf_costs(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calculate_leaf_costs(tree.leafs())


class BrokenBondsScorer(Scorer):
    """Class for scoring nodes and reaction trees based on the breaking of atom bonds

    The score is a summation of the depths in the tree where the focussed bonds
    are found to break in the reaction. If a focussed bond is found to be unbroken
    in the entire tree, the total length of the tree will be added to the score.
    """

    scorer_name = "broken bonds"

    def __init__(self, config: Configuration) -> None:
        super().__init__(config)
        self._break_bonds = [tuple(bond) for bond in config.search.break_bonds]
        self._break_bonds_operator = config.search.break_bonds_operator.lower()
        self._reverse_order = False

    def __repr__(self) -> str:
        return self.scorer_name

    def _calculate_broken_bonds_score(
        self, reactions: Sequence[RetroReaction], depths: Sequence[int]
    ) -> float:

        if not reactions:
            return 0

        max_score = len(set(depths)) * len(self._break_bonds)
        broken_focussed_bonds = []
        scores = []

        # The score should be 0 when no transformations were reported
        if not depths:
            return 0

        for reaction, depth in zip(reactions, depths):
            broken_bonds = BrokenBonds(self._break_bonds)(reaction)
            broken_untracked_bonds = [
                bond for bond in broken_bonds if bond not in broken_focussed_bonds
            ]
            if len(broken_untracked_bonds) > 0:
                broken_focussed_bonds += broken_untracked_bonds
                scores.append(1 - (depth / max_score))
                if self._break_bonds_operator == "or":
                    break

        if self._break_bonds_operator != "or" or (
            self._break_bonds_operator == "or" and len(broken_focussed_bonds) == 0
        ):
            # type: ignore
            if unbroken_bonds := set(self._break_bonds) - set(broken_focussed_bonds):
                scores.append(1 - (len(set(depths)) * len(unbroken_bonds)) / max_score)
        return sum(scores) / len(scores)

    def _score_node(self, node: MctsNode) -> float:
        reactions = node.actions_to()
        depths = []
        for reaction in reactions:
            depth = (
                max(
                    reactant.transform
                    for reactant in reaction.reactants[reaction.index]
                )
                - 1
            )
            depths.append(depth)
        return self._calculate_broken_bonds_score(reactions, depths)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        reactions = [
            reaction.to_smiles_based_retroreaction() for reaction in tree.reactions()
        ]
        depths = [tree.depth(reaction) // 2 for reaction in tree.reactions()]
        return self._calculate_broken_bonds_score(reactions, depths)


class RouteSimilarityScorer(Scorer):
    """
    Class for scoring based on an LSTM model for computing Tree Edit Distance to
    a set of reference routes.

    :param config: the configuration of the tree search
    :param routes_path: the filename of a JSON file with reference routes
    :param model_path: the filename of a checkpoint file with the LSTM model
    :param scaler_params: the parameter settings of the scaler
    :param agg_func: the name of numpy function used to aggregate the distances
                     to the reference routes
    :param similarity: if True, will compute similarity score else distance scores
    """

    scorer_name = "route similarity"

    def __init__(
        self,
        config: Configuration,
        routes_path: str,
        model_path: str,
        scaler_params: Optional[StrDict] = None,
        agg_func: str = "min",
        similarity: bool = False,
    ) -> None:
        if not SUPPORT_DISTANCES:
            raise ValueError(
                "Distance calculations are not supported by this installation."
                " Please install aizynthfinder with extras dependencies."
            )

        # Default scaler from benchmarking
        if scaler_params is None:
            scaler_params = {
                "name": "squash",
                "slope": 0.5,
                "xoffset": 10,
                "yoffset": 0,
            }

        super().__init__(config, scaler_params)
        self._reverse_order = similarity
        self.calculator = route_distances_calculator("lstm", model_path=model_path)

        # Scalar needs to be applied before taking into account of `similarity` parameter,
        # hence we are reassigning the global scaler and then disabling it.
        self._local_scaler = self._scaler
        self._scaler = None

        if not hasattr(np, agg_func):
            raise ValueError(f"Cannot identify aggregate function {agg_func} in numpy")
        self.agg_func = getattr(np, agg_func)
        self.similarity = similarity

        try:
            with open(routes_path or "") as file:
                self.routes = json.load(file)
        except FileNotFoundError:
            logger().info(
                f"Could not load reference routes from {routes_path}. Assuming they will be set later"
            )
            self.routes = []
        self.n_routes = len(self.routes)

    def _score_node(self, node: MctsNode) -> float:
        # We don't have any short-cut to score a node,
        # so we need to convert it to a reaction tree
        return self._score_reaction_tree(node.to_reaction_tree())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        if not self.routes:
            return 0.0 if self.similarity else 1.0

        dist_matrix = self.calculator(self.routes + [tree.to_dict()])
        distances = dist_matrix[self.n_routes, : self.n_routes]
        score = self.agg_func(distances)
        norm_score = self._local_scaler(score)
        if self.similarity:
            return 1.0 - float(norm_score)
        return float(norm_score)


class DeltaSyntheticComplexityScorer(Scorer):
    """
    Class for scoring nodes based on the delta-synthetic-complexity of the node
    and its parent 'horizon' steps up in the tree.

    :param config: the configuration the tree search
    :param sc_score_model: the path to the SCScore model
    :param scaler_params: the parameter settings of the scaler, defaults to max-min between -1.5 and 4
    :param horizon: the number of steps backwards to look for parent molecule
    """

    scorer_name: str = "delta-SC score"

    def __init__(
        self,
        config: Configuration,
        sc_score_model: str,
        scaler_params: Optional[StrDict] = None,
        horizon: int = 3,
    ) -> None:
        if scaler_params is None:
            scaler_params = {
                "name": "min_max",
                "min_val": -1.5,
                "max_val": 4,
                "reverse": False,
            }
        super().__init__(config, scaler_params)
        ## This is necessary because config should not be optional for this scorer
        self._config: Configuration = config

        self.horizon = horizon
        self._model = SCScore(sc_score_model)

    def sc_deltas(self, mols: _Molecules, parents: _Molecules) -> Sequence[float]:
        """
        Calculate delta SC-score among the list of mol-parent pairs.

        :params mols: the leaves of the tree
        :param parents: the parent of the leaves
        :returns: the pair-wise difference in SCScore
        """
        delta_sc_scores = []
        for mol, parent in zip(mols, parents):
            mol.sanitize()
            parent.sanitize()

            sc_score_mol = self._model(mol.rd_mol)
            sc_score_parent = self._model(parent.rd_mol)

            delta_sc_score = sc_score_parent - sc_score_mol
            delta_sc_scores.append(delta_sc_score)
        return delta_sc_scores

    def _get_parent_from_reaction_tree(
        self, tree: ReactionTree, mol: UniqueMolecule, horizon: int
    ) -> UniqueMolecule:
        horizon = min(horizon, tree.depth(mol) // 2)

        parent = mol
        # Step up in the tree, first a reaction, then a molecule
        for _ in range(horizon):
            parent = tree.parent_molecule(parent)

        return parent

    def _get_parent_from_tree_molecule(
        self, mol: TreeMolecule, horizon: int
    ) -> TreeMolecule:
        horizon = min(horizon, mol.transform)

        parent = mol
        # Step up in the tree via parent nodes
        for _ in range(horizon):
            grandparent = parent.parent
            if not grandparent:
                break
            parent = grandparent

        return parent

    def _score_node(self, node: MctsNode) -> float:
        expandable_mols = node.state.expandable_mols
        if not expandable_mols:
            expandable_mols = list(node.state.mols)

        parent_mols = [
            self._get_parent_from_tree_molecule(mol, self.horizon)
            for mol in expandable_mols
        ]
        delta_sc_scores = self.sc_deltas(expandable_mols, parent_mols)
        return min(delta_sc_scores)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        expandable_mols = [node for node in tree.leafs() if not tree.in_stock(node)]
        if not expandable_mols:
            expandable_mols = list(tree.leafs())

        parent_mols = [
            self._get_parent_from_reaction_tree(tree, node, self.horizon)
            for node in expandable_mols
        ]
        delta_sc_scores = self.sc_deltas(expandable_mols, parent_mols)
        return min(delta_sc_scores)


class CombinedScorer(Scorer):
    """Class for scoring nodes and reaction trees by combining weighted scores from a list of scorers

    If no weights are provided as input, the scorer provides default weights that are
    equal for all input scorers.

    The CombinedScorer cannot be instantiated from the config file as it requires the
    names of the scorers to combine as input.
    """

    def __init__(
        self,
        config: Configuration,
        scorers: Sequence[str],
        weights: Optional[Sequence[float]] = None,
    ) -> None:
        super().__init__(config)
        self._scorers = [config.scorers[scorer] for scorer in scorers]
        self._weights = weights if weights else [1 / len(scorers)] * len(scorers)

    def __repr__(self) -> str:
        return " + ".join([repr(scorer) for scorer in self._scorers])

    def _combine_score(self, scores: Sequence[float]) -> float:
        return sum(
            score * weight for score, weight in zip(scores, self._weights)
        ) / sum(self._weights)

    def _score_node(self, node: MctsNode) -> float:
        scores = [scorer(node) for scorer in self._scorers]
        return self._combine_score(scores)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        scores = [scorer(tree) for scorer in self._scorers]
        return self._combine_score(scores)
