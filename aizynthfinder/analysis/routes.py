""" Module containing classes to store and manipulate collections of synthetic routes.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from route_distances.clustering import ClusteringHelper
from route_distances.route_distances import route_distances_calculator

from aizynthfinder.analysis.utils import (
    CombinedReactionTrees,
    RouteSelectionArguments,
)
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsSearchTree, MctsNode
from aizynthfinder.analysis import TreeAnalysis

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        PilImage,
        Optional,
        Any,
        Dict,
        Sequence,
    )
    from aizynthfinder.context.scoring import Scorer


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
    :ivar nodes: the top-ranked MCTS-like nodes
    :ivar scores: initial scores of top-ranked nodes or routes
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
    def from_analysis(
        cls, analysis: TreeAnalysis, selection: RouteSelectionArguments = None
    ) -> "RouteCollection":
        """
        Create a collection from a tree analysis.

        :param analysis: the tree analysis to use
        :param selection: selection criteria for the routes
        :return: the created collection
        """
        items, scores = analysis.sort(selection)
        all_scores = [{repr(analysis.scorer): score} for score in scores]
        kwargs = {"scores": scores, "all_scores": all_scores}
        if isinstance(analysis.search_tree, MctsSearchTree):
            kwargs["nodes"] = items
            reaction_trees = [
                from_node.to_reaction_tree()
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
        cache_key = kwargs.get("model_path", "") if model == "lstm" else content
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
