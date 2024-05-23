""" Module containing classes to perform analysis of the tree search results.
"""

from __future__ import annotations

import operator
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from paretoset import paretorank, paretoset

from aizynthfinder.analysis.utils import RouteSelectionArguments
from aizynthfinder.chem import FixedRetroReaction, hash_reactions
from aizynthfinder.context.scoring import StateScorer
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.andor_trees import AndOrSearchTreeBase
from aizynthfinder.search.mcts import MctsNode, MctsSearchTree

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.context.scoring import Scorer
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        Iterable,
        List,
        Optional,
        Sequence,
        StrDict,
        Tuple,
        Union,
    )

    _Solution = Union[MctsNode, ReactionTree]
    _AnyListOfSolutions = Union[Sequence[MctsNode], Sequence[ReactionTree]]
    _AnyListOfReactions = Sequence[
        Union[Iterable[RetroReaction], Iterable[FixedRetroReaction]]
    ]


class TreeAnalysis:
    """
    Class that encapsulate various analysis that can be
    performed on a search tree.

    :ivar scorers: the objects used to score the nodes
    :ivar search_tree: the search tree

    :param search_tree: the search tree to do the analysis on
    :param scorer: the object used to score the nodes, defaults to StateScorer
    """

    def __init__(
        self,
        search_tree: Union[MctsSearchTree, AndOrSearchTreeBase],
        scorer: Optional[Union[Scorer, List[Scorer]]] = None,
    ) -> None:
        self.search_tree = search_tree
        scorer = scorer or StateScorer(search_tree.config)
        if isinstance(scorer, list):
            self.scorers: List[Scorer] = scorer
        else:
            self.scorers = [scorer]
        self._direction = "max"  # current implementation assumes maximisation
        self._single_objective = len(self.scorers) == 1

    def best(self) -> _Solution:
        """
        Returns the route or MCTS-like node with the highest score.
        If several routes have the same score, it will return the first

        :return: the top scoring node or route
        :raises ValueError: if this is a multi-objective analysis, a single best solution cannot be return
        """
        if not self._single_objective:
            raise ValueError("Cannot return best item for multi-objective analysis")

        if isinstance(self.search_tree, MctsSearchTree):
            nodes = self._all_nodes()
            sorted_nodes, _, _ = self.scorers[0].sort(nodes)
            return sorted_nodes[0]

        sorted_routes, _, _ = self.scorers[0].sort(self.search_tree.routes())
        return sorted_routes[0]

    def pareto_front(self) -> Tuple[_Solution, ...]:
        """
        Returns the Routes or MCTS-like Nodes that are the Pareto Front in the Multi-Objective search.

        :returns: the pareto front solutions.
        :raises ValueError: if this is a single-objective analysis, there is no Pareto-front
        """
        if self._single_objective:
            raise ValueError("Cannot return Pareto front for single-objective analysis")

        if isinstance(self.search_tree, MctsSearchTree):
            solutions = self.search_tree.nodes()
        else:
            solutions = self.search_tree.routes()  # type: ignore

        scores_arr = np.array(
            [[scorer(solution) for scorer in self.scorers] for solution in solutions]
        )
        direction_arr = np.repeat(self._direction, len(self.scorers))
        pareto_mask = paretoset(scores_arr, sense=direction_arr, distinct=False)
        pareto_idxs = np.arange(len(solutions))[pareto_mask]
        pareto_front: Sequence[_Solution] = [solutions[idx] for idx in pareto_idxs]

        return tuple(pareto_front)

    def sort(
        self, selection: Optional[RouteSelectionArguments] = None
    ) -> Tuple[_AnyListOfSolutions, Sequence[Dict[str, float]]]:
        """
        Sort and select the nodes or routes in the search tree.

        The score for each solution is a dictionary of scores,
        one for each of the objectives.

        :param selection: selection criteria for the routes
        :return: the sorted and selected items
        :return: the scores of the sorted items
        """
        selection = selection or RouteSelectionArguments()

        if not self._single_objective:
            return self._pareto_rank_sort(selection)

        scorer = self.scorers[0]
        if isinstance(self.search_tree, MctsSearchTree):
            nodes = self._all_nodes()
            sorted_items, sorted_scores, _ = scorer.sort(nodes)
            actions = [node.actions_to() for node in sorted_items]

        else:
            sorted_items, sorted_scores, _ = scorer.sort(self.search_tree.routes())
            actions = [route.reactions() for route in sorted_items]

        scores = [{repr(scorer): score} for score in sorted_scores]

        return self._collect_top_items(
            sorted_items, scores, actions, sorted_scores, selection, True
        )

    def tree_statistics(self) -> StrDict:
        """
        Returns statistics of the tree

        Currently it returns the number of nodes, the maximum number of transforms,
        maximum number of children, top score, if the top score route is solved,
        the number of molecule in the top score node, and information on pre-cursors

        :return: the statistics
        """
        if isinstance(self.search_tree, MctsSearchTree):
            return self._tree_statistics_mcts()
        return self._tree_statistics_andor()

    def _all_nodes(self) -> Sequence[MctsNode]:
        assert isinstance(self.search_tree, MctsSearchTree)
        # This is to keep backwards compatibility, this should be investigate further
        if repr(self.scorers[0]) == "state score":
            return list(self.search_tree.graph())
        return [node for node in self.search_tree.graph() if not node.children]

    def _pareto_rank_sort(
        self,
        selection: RouteSelectionArguments,
    ) -> Tuple[_AnyListOfSolutions, Sequence[Dict[str, float]]]:
        if isinstance(self.search_tree, MctsSearchTree):
            solutions = self._all_nodes()
        else:
            solutions = self.search_tree.routes()  # type: ignore

        scores_arr = np.array(
            [[scorer(solution) for scorer in self.scorers] for solution in solutions]
        )
        direction_arr = np.repeat(self._direction, len(self.scorers))
        pareto_ranks = paretorank(scores_arr, sense=direction_arr, distinct=False)

        sortidx = sorted(range(len(pareto_ranks)), key=pareto_ranks.__getitem__)
        sorted_pareto_ranks = sorted(pareto_ranks)
        sorted_scores = [
            {
                repr(scorer): scores_arr[idx, scorer_idx]
                for scorer_idx, scorer in enumerate(self.scorers)
            }
            for idx in sortidx
        ]
        sorted_items = [solutions[idx] for idx in sortidx]

        if isinstance(self.search_tree, MctsSearchTree):
            actions = [node.actions_to() for node in sorted_items]
        else:
            actions = [route.reactions() for route in sorted_items]  # type: ignore

        return self._collect_top_items(
            sorted_items, sorted_scores, actions, sorted_pareto_ranks, selection, False
        )

    def _top_nodes(self) -> Tuple[_Solution, ...]:
        if self._single_objective:
            return (self.best(),)
        return self.pareto_front()

    def _tree_statistics_andor(self) -> StrDict:
        assert isinstance(self.search_tree, AndOrSearchTreeBase)
        top_routes = self._top_nodes()
        mols_in_stock = self._top_ranked_join(
            ", ".join(mol.smiles for mol in route.leafs() if route.in_stock(mol))  # type: ignore
            for route in top_routes
        )
        mols_not_in_stock = self._top_ranked_join(
            ", ".join(mol.smiles for mol in route.leafs() if not route.in_stock(mol))  # type: ignore
            for route in top_routes
        )
        all_routes = self.search_tree.routes()
        policy_used_counts = self._policy_used_statistics(
            [reaction for route in all_routes for reaction in route.reactions()]
        )
        availability = self._top_ranked_join(
            ";".join(
                self.search_tree.config.stock.availability_string(mol)
                for mol in route.leafs()  # type: ignore
            )
            for route in top_routes
        )
        number_of_precursors_in_stock = self._top_ranked_join(
            sum(route.in_stock(leaf) for leaf in route.leafs()) for route in top_routes  # type: ignore
        )

        # For multi-objective it becomes to messy to have all the scores in this info,
        # they are nevertheless included with the routes
        if self._single_objective:
            top_score = self.scorers[0](top_routes[0])
        else:
            top_score = None

        return {
            "number_of_nodes": len(self.search_tree.mol_nodes),
            "max_transforms": max(
                node.prop["mol"].transform for node in self.search_tree.mol_nodes
            ),
            "max_children": max(
                len(node.children) for node in self.search_tree.mol_nodes
            ),
            "number_of_routes": len(all_routes),
            "number_of_solved_routes": sum(route.is_solved for route in all_routes),
            "top_score": top_score,
            "is_solved": self._top_ranked_join(route.is_solved for route in top_routes),
            "number_of_steps": self._top_ranked_join(
                len(list(route.reactions())) for route in top_routes  # type: ignore
            ),
            "number_of_precursors": self._top_ranked_join(
                len(list(route.leafs())) for route in top_routes  # type: ignore
            ),
            "number_of_precursors_in_stock": number_of_precursors_in_stock,
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "precursors_availability": availability,
            "policy_used_counts": policy_used_counts,
            "profiling": getattr(self.search_tree, "profiling", {}),
        }

    def _tree_statistics_mcts(self) -> StrDict:
        assert isinstance(self.search_tree, MctsSearchTree)
        top_nodes = self._top_nodes()
        assert isinstance(top_nodes[0], MctsNode)
        top_states = [node.state for node in top_nodes]  # type: ignore
        nodes = list(self.search_tree.graph())
        mols_in_stock = self._top_ranked_join(
            ", ".join(
                mol.smiles
                for mol, instock in zip(state.mols, state.in_stock_list)
                if instock
            )
            for state in top_states
        )
        mols_not_in_stock = self._top_ranked_join(
            ", ".join(
                mol.smiles
                for mol, instock in zip(state.mols, state.in_stock_list)
                if not instock
            )
            for state in top_states
        )

        policy_used_counts = self._policy_used_statistics(
            [node[child]["action"] for node in nodes for child in node.children]
        )

        # For multi-objective it becomes to messy to have all the scores in this info,
        # they are nevertheless included with the routes
        if self._single_objective:
            top_score = self.scorers[0](top_nodes[0])
        else:
            top_score = None

        return {
            "number_of_nodes": len(nodes),
            "max_transforms": max(node.state.max_transforms for node in nodes),
            "max_children": max(len(node.children) for node in nodes),
            "number_of_routes": sum(1 for node in nodes if not node.children),
            "number_of_solved_routes": sum(
                1 for node in nodes if not node.children and node.state.is_solved
            ),
            "top_score": top_score,
            "is_solved": self._top_ranked_join(state.is_solved for state in top_states),
            "number_of_steps": self._top_ranked_join(
                state.max_transforms for state in top_states
            ),
            "number_of_precursors": self._top_ranked_join(
                len(state.mols) for state in top_states
            ),
            "number_of_precursors_in_stock": self._top_ranked_join(
                sum(state.in_stock_list) for state in top_states
            ),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "precursors_availability": self._top_ranked_join(
                ";".join(state.stock_availability) for state in top_states
            ),
            "policy_used_counts": policy_used_counts,
            "profiling": getattr(self.search_tree, "profiling", {}),
        }

    @staticmethod
    def _collect_top_items(
        items: _AnyListOfSolutions,
        scores: Sequence[Dict[str, float]],
        reactions: _AnyListOfReactions,
        ranks: Sequence[Union[int, float]],
        selection: RouteSelectionArguments,
        min_comp: bool,
    ) -> Tuple[_AnyListOfSolutions, Sequence[Dict[str, float]]]:
        """
        Finds and returns the top-ranked super-nodes or reaction trees and
        their scores.

        The `selection` argument is used to control the number of items that are returned.
        If at least one of the items is solved, and `selection` specifies that all items
        should be returned, all solved items are returned.
        Otherwise, at least a minimum number of items are returned, but there could be more
        than this if they have the same rank. However not more than a maximum number of
        items are returned.

        Duplicated routes are ignored.

        The items are compared based on rank. For single-ojective analysis, this is simply
        the score and higher is better. For multi-objective analysis, the rank is the Pareto
        rank and smaller i better.

        :param items: the super-nodes or reaction trees to select from, sorted by rank
        :param scores: the computed scores of the nodes
        :param reaction: the reactions for each item, used for checking duplicity
        :param ranks: the rank of each item, used for determining selection
        :param selection: min and max solutions, or indication to return all
        :param min_comp: True if higher rank is better, else False
        """
        if len(items) <= selection.nmin:
            return items, scores

        max_return, min_return = selection.nmax, selection.nmin
        if selection.return_all:
            nsolved = sum(int(item.is_solved) for item in items)
            if nsolved:
                max_return = nsolved
                min_return = nsolved

        seen_hashes = set()
        best_items: List[Any] = []
        best_scores = []
        if min_comp:
            last_rank = 1e16
            comp_op = operator.lt
        else:
            last_rank = 0
            comp_op = operator.gt
        for rank, score, item, actions in zip(ranks, scores, items, reactions):
            if len(best_items) >= min_return and comp_op(rank, last_rank):
                break
            route_hash = hash_reactions(actions)

            if route_hash in seen_hashes:
                continue
            seen_hashes.add(route_hash)
            best_items.append(item)
            best_scores.append(score)
            last_rank = rank

            if max_return and len(best_items) == max_return:
                break

        return best_items, best_scores

    @staticmethod
    def _policy_used_statistics(
        reactions: Iterable[Union[RetroReaction, FixedRetroReaction]]
    ) -> StrDict:
        policy_used_counts: StrDict = defaultdict(lambda: 0)
        for reaction in reactions:
            policy_used = reaction.metadata.get("policy_name")
            if policy_used:
                policy_used_counts[policy_used] += 1
        return dict(policy_used_counts)

    @staticmethod
    def _top_ranked_join(items: Iterable[Any]) -> Union[str, Any]:
        items = list(items)
        if len(items) == 1:
            return items[0]
        return "|".join(f"{item}" for item in items)
