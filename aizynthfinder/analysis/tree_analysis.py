""" Module containing classes to perform analysis of the tree search results.
"""
from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING

from aizynthfinder.chem import (
    FixedRetroReaction,
    hash_reactions,
)
from aizynthfinder.context.scoring import (
    StateScorer,
)
from aizynthfinder.analysis.utils import RouteSelectionArguments
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsSearchTree, MctsNode
from aizynthfinder.search.andor_trees import AndOrSearchTreeBase

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        Union,
        Tuple,
        Any,
        Iterable,
        Sequence,
        List,
    )
    from aizynthfinder.context.scoring import Scorer
    from aizynthfinder.chem import RetroReaction


class TreeAnalysis:
    """
    Class that encapsulate various analysis that can be
    performed on a search tree.

    :ivar scorer: the object used to score the nodes
    :ivar search_tree: the search tree

    :param search_tree: the search tree to do the analysis on
    :param scorer: the object used to score the nodes, defaults to StateScorer
    """

    def __init__(
        self,
        search_tree: Union[MctsSearchTree, AndOrSearchTreeBase],
        scorer: Scorer = None,
    ) -> None:
        self.search_tree = search_tree
        if scorer is None:
            self.scorer: Scorer = StateScorer(search_tree.config)
        else:
            self.scorer = scorer

    def best(self) -> Union[MctsNode, ReactionTree]:
        """
        Returns the route or MCTS-like node with the highest score.
        If several routes have the same score, it will return the first

        :return: the top scoring node or route
        """
        if isinstance(self.search_tree, MctsSearchTree):
            nodes = self._all_nodes()
            sorted_nodes, _, _ = self.scorer.sort(nodes)
            return sorted_nodes[0]

        sorted_routes, _, _ = self.scorer.sort(self.search_tree.routes())
        return sorted_routes[0]

    def sort(
        self, selection: RouteSelectionArguments = None
    ) -> Tuple[Union[Sequence[MctsNode], Sequence[ReactionTree]], Sequence[float]]:
        """
        Sort and select the nodes or routes in the search tree.

        :param selection: selection criteria for the routes
        :return: the items
        :return: the score
        """
        selection = selection or RouteSelectionArguments()

        if isinstance(self.search_tree, MctsSearchTree):
            nodes = self._all_nodes()
            sorted_items, sorted_scores, _ = self.scorer.sort(nodes)
            actions = [node.actions_to() for node in sorted_items]

        else:
            sorted_items, sorted_scores, _ = self.scorer.sort(self.search_tree.routes())
            actions = [route.reactions() for route in sorted_items]

        return self._collect_top_items(sorted_items, sorted_scores, actions, selection)

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
        if repr(self.scorer) == "state score":
            return list(self.search_tree.graph())
        return [node for node in self.search_tree.graph() if not node.children]

    def _tree_statistics_andor(self) -> StrDict:
        assert isinstance(self.search_tree, AndOrSearchTreeBase)
        top_route = self.best()
        assert isinstance(top_route, ReactionTree)
        mols_in_stock = ", ".join(
            mol.smiles for mol in top_route.leafs() if top_route.in_stock(mol)
        )
        mols_not_in_stock = ", ".join(
            mol.smiles for mol in top_route.leafs() if not top_route.in_stock(mol)
        )
        all_routes = self.search_tree.routes()
        policy_used_counts = self._policy_used_statistics(
            [reaction for route in all_routes for reaction in route.reactions()]
        )
        availability = ";".join(
            self.search_tree.config.stock.availability_string(mol)
            for mol in top_route.leafs()
        )

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
            "top_score": self.scorer(top_route),
            "is_solved": top_route.is_solved,
            "number_of_steps": len(list(top_route.reactions())),
            "number_of_precursors": len(list(top_route.leafs())),
            "number_of_precursors_in_stock": sum(
                top_route.in_stock(leaf) for leaf in top_route.leafs()
            ),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "precursors_availability": availability,
            "policy_used_counts": policy_used_counts,
        }

    def _tree_statistics_mcts(self) -> StrDict:
        assert isinstance(self.search_tree, MctsSearchTree)
        top_node = self.best()
        assert isinstance(top_node, MctsNode)
        top_state = top_node.state
        nodes = list(self.search_tree.graph())
        mols_in_stock = ", ".join(
            mol.smiles
            for mol, instock in zip(top_state.mols, top_state.in_stock_list)
            if instock
        )
        mols_not_in_stock = ", ".join(
            mol.smiles
            for mol, instock in zip(top_state.mols, top_state.in_stock_list)
            if not instock
        )

        policy_used_counts = self._policy_used_statistics(
            [node[child]["action"] for node in nodes for child in node.children]
        )

        return {
            "number_of_nodes": len(nodes),
            "max_transforms": max(node.state.max_transforms for node in nodes),
            "max_children": max(len(node.children) for node in nodes),
            "number_of_routes": sum(1 for node in nodes if not node.children),
            "number_of_solved_routes": sum(
                1 for node in nodes if not node.children and node.state.is_solved
            ),
            "top_score": self.scorer(top_node),
            "is_solved": top_state.is_solved,
            "number_of_steps": top_state.max_transforms,
            "number_of_precursors": len(top_state.mols),
            "number_of_precursors_in_stock": sum(top_state.in_stock_list),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "precursors_availability": ";".join(top_state.stock_availability),
            "policy_used_counts": policy_used_counts,
        }

    @staticmethod
    def _collect_top_items(
        items: Union[Sequence[MctsNode], Sequence[ReactionTree]],
        scores: Sequence[float],
        reactions: Sequence[
            Union[Iterable[RetroReaction], Iterable[FixedRetroReaction]]
        ],
        selection,
    ) -> Tuple[Union[Sequence[MctsNode], Sequence[ReactionTree]], Sequence[float]]:
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
        last_score = 1e16
        for score, item, actions in zip(scores, items, reactions):
            if len(best_items) >= min_return and score < last_score:
                break
            route_hash = hash_reactions(actions)

            if route_hash in seen_hashes:
                continue
            seen_hashes.add(route_hash)
            best_items.append(item)
            best_scores.append(score)
            last_score = score

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
