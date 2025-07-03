""" Module containing classes used to score the reaction routes.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from rxnutils.chem.features.sc_score import SCScore
from rxnutils.routes.scoring import DeepsetModelClient, deepset_route_score

try:
    from route_distances.route_distances import route_distances_calculator
except ImportError:
    SUPPORT_DISTANCES = False
else:
    SUPPORT_DISTANCES = True

from aizynthfinder.context.scoring.scorers_base import Scorer, make_rxnutils_route

# pylint: disable=W0611
from aizynthfinder.context.scoring.scorers_mols import (
    DeltaSyntheticComplexityScorer,
    FractionInSourceStockScorer,
    FractionInStockScorer,
    FractionOfIntermediatesInStockScorer,
    NumberOfPrecursorsInStockScorer,
    NumberOfPrecursorsScorer,
    PriceSumScorer,
    StockAvailabilityScorer,
)
from aizynthfinder.context.scoring.scorers_reactions import (
    AverageTemplateOccurrenceScorer,
    MaxTransformScorer,
    NumberOfReactionsScorer,
    ReactionClassMembershipScorer,
    ReactionClassRankScorer,
)
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode
from aizynthfinder.utils.bonds import BrokenBonds
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Optional, Sequence, StrDict, TypeVar

    _Scoreable = TypeVar("_Scoreable", MctsNode, ReactionTree)


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
        combine_strategy: str = "mean-arithmetic",
        short_name: Optional[str] = None,
    ) -> None:
        super().__init__(config)
        self._scorers = [config.scorers[scorer] for scorer in scorers]
        self._weights = weights if weights else [1 / len(scorers)] * len(scorers)

        if combine_strategy not in ["mean-arithmetic", "mean-geometric", "product"]:
            raise ValueError(
                f"'combine_strategy' should be one of ['mean-arithmetic', "
                f"'mean-geometric', 'product'], got {combine_strategy}"
            )
        self._combine_strategy = combine_strategy
        self._short_name = short_name

    def __repr__(self) -> str:
        if self._short_name:
            return self._short_name
        scorer_name = " + ".join([repr(scorer) for scorer in self._scorers])
        if self._combine_strategy == "mean-arithmetic":
            return scorer_name
        return f"{scorer_name} ({self._combine_strategy})"

    def _combine_score(self, scores: Sequence[float]) -> float:
        if self._combine_strategy == "mean-arithmetic":
            return sum(
                score * weight for score, weight in zip(scores, self._weights)
            ) / sum(self._weights)

        product_score = 1.0
        for score in scores:
            product_score *= score

        if self._combine_strategy == "mean-geometric":
            return product_score ** (1 / len(scores))
        return product_score

    def _score_node(self, node: MctsNode) -> float:
        scores = [scorer(node) for scorer in self._scorers]
        return self._combine_score(scores)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        scores = [scorer(tree) for scorer in self._scorers]
        return self._combine_score(scores)


class StateScorer(Scorer):
    """Class for scoring nodes based on the state score"""

    scorer_name = "state score"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        # This is necessary because config should not be optional for this scorer
        self._config: Configuration = config
        self._transform_scorer = MaxTransformScorer(
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


class DeepSetScorer(Scorer):
    """Class for scoring nodes and reaction trees by a deep learning model,
    augmented by expert knowledge.

    Wrapper for `deepset_route_score` from `rxnutils` with addition
    of the expert knowledge correction for route length.

    The raw output of the score is capped between 0 and 20.
    0 because sometimes the deepset model return negative values, even though
    it should be confined to >= 0
    20 because anything large than 20 is useless and uninteresting, this
    cap could probably be removed in the future if we retrain it.

    How to interpret the score:
    - ”Good” for scores between 0 and 5
    - ”Plausible” for scores between 5 and 9
    - ”Bad” for scores between 10 and 15

    :param config: the configuration the tree search
    :param deepset_model: the path to the DeepSet model
    :param sc_score_model: the path to the SCScore model
    :param class_ranks_path: the path to a file with reaction class ranks
    :param scaler_params: the parameter settings of the scaler
    """

    scorer_name = "expert-augmented score"

    def __init__(
        self,
        config: Configuration,
        deepset_model: str,
        sc_score_model: str,
        class_ranks_path: str,
        model_score_weight: float = 0.97,
        length_weight: float = -0.43,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._reverse_order = False

        data = pd.read_csv(class_ranks_path, sep=",")
        self._reaction_class_ranks = dict(
            zip(data["reaction_class"], data["rank_score"])
        )

        self._deepset_model = DeepsetModelClient(deepset_model)
        self._sc_scorer = SCScore(sc_score_model)
        self._model_score_weight = model_score_weight
        self._length_weight = length_weight

    def _score_node(self, node: MctsNode) -> float:
        # All rxnutils scorers are based on reaction trees
        return self._score_reaction_tree(node.to_reaction_tree())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        route = make_rxnutils_route(tree)
        pred_distance = deepset_route_score(
            route, self._deepset_model, self._sc_scorer, self._reaction_class_ranks
        )
        raw_score = (
            self._model_score_weight * pred_distance
            + self._length_weight * route.nsteps
        )
        return min(max(0, raw_score), 20)
