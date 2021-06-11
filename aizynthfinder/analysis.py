""" Module containing classes to perform analysis of the tree search results.
"""
from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from deprecated import deprecated
from route_distances.clustering import ClusteringHelper
from route_distances.route_distances import route_distances_calculator

from aizynthfinder.chem import (
    FixedRetroReaction,
    hash_reactions,
)
from aizynthfinder.utils.analysis_helpers import (
    CombinedReactionTrees,
    ReactionTreeFromMcts,
)
from aizynthfinder.context.scoring import (
    StateScorer,
)
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.mcts.mcts import SearchTree as MctsSearchTree
from aizynthfinder.mcts.node import Node as MctsNode
from aizynthfinder.utils.trees import AndOrSearchTreeBase

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        PilImage,
        Optional,
        Union,
        Tuple,
        Any,
        Iterable,
        Dict,
        Sequence,
        List,
    )
    from aizynthfinder.context.scoring import Scorer
    from aizynthfinder.chem import RetroReaction


@deprecated(reason="use a factory method directly", version="2.5.0")
def rt_from_analysis(
    analysis: TreeAnalysis, from_node: MctsNode = None
) -> "ReactionTree":
    """
    Create a reaction from a tree analysis.

    The single route can be from a given node in the search tree if the `from_node`
    argument is given. If it is not given, the top scoring node is used.

    :param analysis: the analysis to base the reaction tree on
    :param from_node: the end node of the route, defaults to None
    :returns: the reaction tree
    """
    if not isinstance(analysis.search_tree, MctsSearchTree):
        raise TreeAnalysisException(
            "Cannot create reaction tree because analysis not on a MCTS tree"
        )

    if not from_node:
        from_node = analysis.best_node()

    actions, nodes = analysis.search_tree.route_to_node(from_node)
    return ReactionTreeFromMcts(actions=actions, nodes=nodes).tree


ReactionTree.from_analysis = rt_from_analysis  # type: ignore


class TreeAnalysisException(Exception):
    """Exception raised by TreeAnalysis class"""


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
        Returns the route with the highest score.
        If several routes have the same score, it will return the first

        :return: the top scoring node or route
        """
        if isinstance(self.search_tree, MctsSearchTree):
            return self.best_node()

        sorted_routes, _, _ = self.scorer.sort(self.search_tree.routes())
        return sorted_routes[0]

    def best_node(self) -> MctsNode:
        """
        Returns the node with the highest score.
        If several nodes have the same score, it will return the first

        :return: the top scoring node
        :raises TreeAnalysisException: when the search is not a MCTS tree
        """
        if not isinstance(self.search_tree, MctsSearchTree):
            raise TreeAnalysisException(
                "Cannot compute best node because analysis not on a MCTS tree"
            )

        nodes = self._all_nodes()
        sorted_nodes, _, _ = self.scorer.sort(nodes)
        return sorted_nodes[0]

    def sort(
        self, min_return: int = 5
    ) -> Tuple[Union[Sequence[MctsNode], Sequence[ReactionTree]], Sequence[float]]:
        """
        Sort and select the nodes or routes in the search tree.
        The algorithm filter away identical routes and returns at minimum the number specified.
        If multiple alternative routes have the same score as the n'th route, they will be included and returned.


        :param min_return: the minium number of routes to return, defaults to 5
        :return: the items
        :return: the score
        """
        if isinstance(self.search_tree, MctsSearchTree):
            return self.sort_nodes(min_return)

        sorted_routes, sorted_scores, _ = self.scorer.sort(self.search_tree.routes())
        actions = [route.reactions() for route in sorted_routes]
        return self._collect_top_items(
            sorted_routes, sorted_scores, actions, min_return
        )

    def sort_nodes(
        self, min_return: int = 5, max_return: int = 25
    ) -> Tuple[Sequence[MctsNode], Sequence[float]]:
        """
        Sort and select the nodes, so that the best scoring routes are returned.
        The algorithm filter away identical routes and returns at minimum the number specified.
        If multiple alternative routes have the same score as the n'th route, they will be included and returned.

        :param min_return: the minimum number of routes to return, defaults to 5
        :param max_return: the maximum number of routes to return
        :return: the nodes
        :return: the score
        :raises TreeAnalysisException: when the search is not a MCTS tree
        """
        if not isinstance(self.search_tree, MctsSearchTree):
            raise TreeAnalysisException(
                "Cannot sort nodes because analysis not on a MCTS tree"
            )

        nodes = self._all_nodes()
        sorted_nodes, sorted_scores, _ = self.scorer.sort(nodes)
        actions = [self.search_tree.route_to_node(node)[0] for node in sorted_nodes]

        return self._collect_top_items(  # type: ignore
            sorted_nodes, sorted_scores, actions, min_return, max_return
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
        if repr(self.scorer) == "state score":
            return list(self.search_tree.graph())
        return [node for node in self.search_tree.graph() if not node.children()]

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
            "policy_used_counts": policy_used_counts,
        }

    def _tree_statistics_mcts(self) -> StrDict:
        assert isinstance(self.search_tree, MctsSearchTree)
        top_node = self.best_node()
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
            [node[child]["action"] for node in nodes for child in node.children()]
        )

        return {
            "number_of_nodes": len(nodes),
            "max_transforms": max(node.state.max_transforms for node in nodes),
            "max_children": max(len(node.children()) for node in nodes),
            "number_of_routes": sum(1 for node in nodes if not node.children()),
            "number_of_solved_routes": sum(
                1 for node in nodes if not node.children() and node.state.is_solved
            ),
            "top_score": self.scorer(top_node),
            "is_solved": top_state.is_solved,
            "number_of_steps": top_state.max_transforms,
            "number_of_precursors": len(top_state.mols),
            "number_of_precursors_in_stock": sum(top_state.in_stock_list),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "policy_used_counts": policy_used_counts,
        }

    @staticmethod
    def _collect_top_items(
        items: Union[Sequence[MctsNode], Sequence[ReactionTree]],
        scores: Sequence[float],
        reactions: Sequence[
            Union[Iterable[RetroReaction], Iterable[FixedRetroReaction]]
        ],
        min_return: int,
        max_return: int = None,
    ) -> Tuple[Union[Sequence[MctsNode], Sequence[ReactionTree]], Sequence[float]]:
        if len(items) <= min_return:
            return items, scores

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


class RouteCollection:
    """
    Holds a collections of reaction routes.

    If can be the top scored nodes, their scores and
    the reaction trees created from them.
    It can also be a cluster of such routes.

    The class has the functionality to compute collective results
    for the different routes such as images.

    Properties of individual route can be obtained with simple indexing.

    .. code-block::

        route0 = collection[0]

    :ivar all_scores: all the computed scores for the routes
    :ivar nodes: the top-ranked nodes
    :ivar scores: initial scores of top-ranked nodes
    :ivar reaction_trees: the reaction trees created from the top-ranked nodes
    :ivar clusters: the created clusters from the collection

    :param reaction_trees: the trees to base the collection on
    """

    def __init__(self, reaction_trees: Sequence[ReactionTree], **kwargs) -> None:
        self._routes: Sequence[StrDict] = [{} for _ in range(len(reaction_trees))]
        self.reaction_trees = reaction_trees
        self._update_route_dict(reaction_trees, "reaction_tree")

        self.nodes = self._unpack_kwarg_with_default("nodes", None, **kwargs)
        self.scores = self._unpack_kwarg_with_default("scores", np.nan, **kwargs)
        self.all_scores = self._unpack_kwarg_with_default("all_scores", dict, **kwargs)

        self._dicts: Optional[Sequence[StrDict]] = self._unpack_kwarg("dicts", **kwargs)
        self._images: Optional[Sequence[PilImage]] = self._unpack_kwarg(
            "images", **kwargs
        )
        self._jsons: Optional[Sequence[str]] = self._unpack_kwarg("jsons", **kwargs)
        self.clusters: Optional[Sequence[RouteCollection]] = self._unpack_kwarg(
            "clusters", **kwargs
        )
        self._distance_matrix: Dict[str, np.ndarray] = {}
        self._combined_reaction_trees: Optional[CombinedReactionTrees] = None

    @classmethod
    def from_analysis(cls, analysis: TreeAnalysis, min_nodes: int) -> "RouteCollection":
        """
        Create a collection from a tree analysis.

        :param analysis: the tree analysis to use
        :param min_nodes: the minimum number of top-ranked nodes to consider
        :return: the created collection
        """
        items, scores = analysis.sort(min_return=min_nodes)
        all_scores = [{repr(analysis.scorer): score} for score in scores]
        kwargs = {"scores": scores, "all_scores": all_scores}
        if isinstance(analysis.search_tree, MctsSearchTree):
            kwargs["nodes"] = items
            reaction_trees = [
                ReactionTreeFromMcts(
                    *analysis.search_tree.route_to_node(from_node)
                ).tree
                for from_node in items
                if isinstance(from_node, MctsNode)
            ]
        else:
            reaction_trees = items  # type: ignore

        return cls(reaction_trees, **kwargs)

    def __getitem__(self, index: int) -> StrDict:
        if index < 0 or index >= len(self):
            raise IndexError("Index out of range")
        return self._routes[index]

    def __len__(self) -> int:
        return len(self.reaction_trees)

    @property
    def dicts(self) -> Sequence[StrDict]:
        """Returns a list of dictionary representation of the routes"""
        if self._dicts is None:
            self._dicts = self.make_dicts()
        return self._dicts

    @property
    def images(self) -> Sequence[PilImage]:
        """Returns a list of pictoral representation of the routes"""
        if self._images is None:
            self._images = self.make_images()
        return self._images

    @property
    def jsons(self) -> Sequence[str]:
        """Returns a list of JSON string representation of the routes"""
        if self._jsons is None:
            self._jsons = self.make_jsons()
        return self._jsons

    def cluster(
        self,
        n_clusters: int,
        max_clusters: int = 5,
        distances_model: str = "ted",
        **kwargs: Any
    ) -> np.ndarray:
        """
        Cluster the route collection into a number of clusters.

        Additional arguments to the distance or clustering algorithm
        can be passed in as key-word arguments.

        When `distances_model` is "lstm", a key-word argument `model_path` needs to be given
        when `distances_model` is "ted", two optional key-word arguments `timeout` and `content`
        can be given.

        If the number of reaction trees are less than 3, no clustering will be performed

        :param n_clusters: the desired number of clusters, if less than 2 triggers optimization
        :param max_clusters: the maximum number of clusters to consider
        :param distances_model: can be ted or lstm and determines how the route distances are computed
        :return: the cluster labels
        """
        if len(self.reaction_trees) < 3:
            return np.asarray([])
        dist_kwargs = {
            "content": kwargs.pop("content", "both"),
            "timeout": kwargs.pop("timeout", None),
            "model_path": kwargs.pop("model_path", None),
        }
        try:
            distances = self.distance_matrix(model=distances_model, **dist_kwargs)
        except ValueError:
            return np.asarray([])

        labels = ClusteringHelper.cluster(
            distances,
            n_clusters,
            max_clusters=max_clusters,
            **kwargs,
        )
        self._make_clusters(labels)
        return labels

    def combined_reaction_trees(self, recreate: bool = False) -> CombinedReactionTrees:
        """
        Return an object that combines all the reaction tree into a single reaction tree graph

        :param recreate: if False will return a cached object if available, defaults to False
        :return: the combined trees
        """
        if not self._combined_reaction_trees or recreate:
            self._combined_reaction_trees = CombinedReactionTrees(self.reaction_trees)
        return self._combined_reaction_trees

    def compute_scores(self, *scorers: Scorer) -> None:
        """
        Compute new scores for all routes in this collection.
        They can then be accessed with the ``all_scores`` attribute.
        """
        if self.nodes[0]:
            list_ = self.nodes
        else:
            list_ = self.reaction_trees

        for scorer in scorers:
            for idx, score in enumerate(scorer(list_)):  # type: ignore
                self.all_scores[idx][repr(scorer)] = score
        self._update_route_dict(self.all_scores, "all_score")

    def dict_with_scores(self) -> Sequence[StrDict]:
        """
        Return the routes as dictionaries with all scores added
        to the root (target) node.

        :return: the routes as dictionaries
        """
        dicts = []
        for dict_, scores in zip(self.dicts, self.all_scores):
            dicts.append(dict(dict_))
            dicts[-1]["scores"] = dict(scores)
        return dicts

    def distance_matrix(
        self, recreate: bool = False, model: str = "ted", **kwargs: Any
    ) -> np.ndarray:
        """
        Compute the distance matrix between each pair of reaction trees

        All key-word arguments are passed along to the `route_distance_calculator`
        function from the `route_distances` package.

        When `model` is "lstm", a key-word argument `model_path` needs to be given
        when `model` is "ted", two optional key-word arguments `timeout` and `content`
        can be given.

        :param recreate: if False, use a cached one if available
        :param model: the type of model to use "ted" or "lstm"
        :return: the square distance matrix
        """
        if model == "lstm" and not kwargs.get("model_path"):
            raise KeyError(
                "Need to provide 'model_path' argument when using LSTM model for computing distances"
            )
        content = kwargs.get("content", "both")
        cache_key = kwargs.get("model_path") if model == "lstm" else content
        if self._distance_matrix.get(cache_key) is not None and not recreate:
            return self._distance_matrix[cache_key]
        calculator = route_distances_calculator(model, **kwargs)
        distances = calculator(self.dicts)
        self._distance_matrix[cache_key] = distances
        return distances

    def make_dicts(self) -> Sequence[StrDict]:
        """Convert all reaction trees to dictionaries"""
        self._dicts = [tree.to_dict() for tree in self.reaction_trees]
        self._update_route_dict(self._dicts, "dict")
        return self._dicts

    def make_images(self) -> Sequence[Optional[PilImage]]:
        """Convert all reaction trees to images"""

        self._images = []
        for tree in self.reaction_trees:
            try:
                img = tree.to_image()
            except ValueError:
                self._images.append(None)
            else:
                self._images.append(img)
        self._update_route_dict(self._images, "image")
        return self._images

    def make_jsons(self) -> Sequence[str]:
        """Convert all reaction trees to JSON strings"""
        self._jsons = [tree.to_json() for tree in self.reaction_trees]
        self._update_route_dict(self._jsons, "json")
        return self._jsons

    def rescore(self, scorer: Scorer) -> None:
        """
        Rescore the routes in the collection, and thereby re-order them.

        This will replace the ``scores`` attribute, and update the ``all_scores``
        attribute with another entry.

        :param scorer: the scorer to use
        """
        if self.nodes[0]:
            self.nodes, self.scores, sortidx = scorer.sort(self.nodes)
            self.reaction_trees = [self.reaction_trees[idx] for idx in sortidx]
        else:
            self.reaction_trees, self.scores, sortidx = scorer.sort(self.reaction_trees)
        self._routes = [self._routes[idx] for idx in sortidx]
        self.all_scores = [self.all_scores[idx] for idx in sortidx]
        if self._dicts:
            self._dicts = [self._dicts[idx] for idx in sortidx]
        if self._images:
            self._images = [self._images[idx] for idx in sortidx]
        if self._jsons:
            self._jsons = [self._jsons[idx] for idx in sortidx]

        for idx, score in enumerate(self.scores):
            self.all_scores[idx][repr(scorer)] = score
        self._update_route_dict(self.all_scores, "all_score")

    def _make_clusters(self, clusters: np.ndarray) -> None:
        n_clusters = max(clusters) + 1
        self.clusters = []
        for cluster in range(n_clusters):
            selection = clusters == cluster
            kwargs = {
                "reaction_trees": self._select_subset(self.reaction_trees, selection),
                "nodes": self._select_subset(self.nodes, selection),
                "scores": self._select_subset(self.scores, selection),
            }
            if self._images:
                kwargs["images"] = self._select_subset(self.images, selection)
            if self._dicts:
                kwargs["dicts"] = self._select_subset(self.dicts, selection)
            if self._jsons:
                kwargs["jsons"] = self._select_subset(self.jsons, selection)

            self.clusters.append(RouteCollection(**kwargs))

    def _unpack_kwarg(self, key: str, **kwargs: Any) -> Optional[Sequence[Any]]:
        if key not in kwargs:
            return None
        arr = kwargs[key]
        self._update_route_dict(arr, key[:-1])
        return arr

    def _unpack_kwarg_with_default(
        self, key: str, default: Any, **kwargs: Any
    ) -> Sequence[Any]:
        arr = self._unpack_kwarg(key, **kwargs)
        if arr is not None:
            return arr
        return [
            default() if callable(default) else default
            for _ in range(len(self.reaction_trees))
        ]

    def _update_route_dict(self, arr: Sequence[Any], key: str) -> None:
        for i, value in enumerate(arr):
            self._routes[i][key] = value

    @staticmethod
    def _select_subset(arr: Sequence[Any], selection: Sequence[bool]) -> Sequence[Any]:
        return [item for sel, item in zip(selection, arr) if sel]
