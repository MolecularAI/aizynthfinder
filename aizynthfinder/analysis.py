""" Module containing classes to perform analysis of the tree search results.
"""
from __future__ import annotations
import json
import time
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
import networkx as nx

from aizynthfinder.chem import (
    Molecule,
    UniqueMolecule,
    FixedRetroReaction,
    hash_reactions,
    none_molecule,
)
from aizynthfinder.utils.image import make_graphviz_image
from aizynthfinder.utils.analysis_helpers import (
    ReactionTreeFromDict,
    ReactionTreeFromMcts,
    CombinedReactionTrees,
)
from aizynthfinder.utils.route_clustering import ReactionTreeWrapper, ClusteringHelper

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        PilImage,
        FrameColors,
        Optional,
        Union,
        Tuple,
        Any,
        Iterable,
        List,
        Dict,
        Sequence,
    )
    from aizynthfinder.context.scoring import Scorer
    from aizynthfinder.mcts.mcts import SearchTree
    from aizynthfinder.mcts.node import Node


class TreeAnalysis:
    """
    Class that encapsulate various analysis that can be
    performed on a search tree.

    :ivar scorer: the object used to score the nodes
    :ivar search_tree: the search tree

    :param search_tree: the search tree to do the analysis on
    :param scorer: the object used to score the nodes, defaults to StateScorer
    """

    def __init__(self, search_tree: SearchTree, scorer: Scorer = None) -> None:
        self.search_tree = search_tree
        if scorer is None:
            # Do import here to avoid circular imports
            from aizynthfinder.context.scoring import StateScorer

            self.scorer: Scorer = StateScorer(search_tree.config)
        else:
            self.scorer = scorer

    def best_node(self) -> Node:
        """
        Returns the node with the highest score.
        If several nodes have the same score, it will return the first

        :return: the top scoring node
        """
        nodes = self._all_nodes()
        sorted_nodes, _, _ = self.scorer.sort(nodes)
        return sorted_nodes[0]

    def sort_nodes(
        self, min_return: int = 5, max_return: int = 25
    ) -> Tuple[Sequence[Node], Sequence[float]]:
        """
        Sort and select the nodes, so that the best scoring routes are returned.
        The algorithm filter away identical routes and returns at minimum the number specified.
        If multiple alternative routes have the same score as the n'th route, they will be included and returned.

        :param min_return: the minimum number of routes to return, defaults to 5
        :param max_return: the maximum number of routes to return
        :return: the nodes
        :return: the score
        """
        nodes = self._all_nodes()
        sorted_nodes, sorted_scores, _ = self.scorer.sort(nodes)

        if len(nodes) <= min_return:
            return sorted_nodes, sorted_scores

        seen_hashes = set()
        best_nodes: List[Node] = []
        best_scores: List[float] = []
        last_score = 1e16
        for score, node in zip(sorted_scores, sorted_nodes):
            if len(best_nodes) >= min_return and score < last_score:
                break
            route_actions, _ = self.search_tree.route_to_node(node)
            route_hash = hash_reactions(route_actions)

            if route_hash in seen_hashes:
                continue
            seen_hashes.add(route_hash)
            best_nodes.append(node)
            best_scores.append(score)
            last_score = score

            if max_return and len(best_nodes) == max_return:
                break

        return best_nodes, best_scores

    def tree_statistics(self) -> StrDict:
        """
        Returns statistics of the tree

        Currently it returns the number of nodes, the maximum number of transforms,
        maximum number of children, top score, if the top score route is solved,
        the number of molecule in the top score node, and information on pre-cursors

        :return: the statistics
        """
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

        policy_used_counts: StrDict = defaultdict(lambda: 0)
        for node in nodes:
            for child in node.children():
                policy_used = node[child]["action"].metadata.get("policy_name")
                if policy_used:
                    policy_used_counts[policy_used] += 1

        return {
            "number_of_nodes": len(nodes),
            "max_transforms": max(node.state.max_transforms for node in nodes),
            "max_children": max(len(node.children()) for node in nodes),
            "number_of_leafs": sum(1 for node in nodes if not node.children()),
            "number_of_solved_leafs": sum(
                1 for node in nodes if not node.children() and node.state.is_solved
            ),
            "top_score": self.scorer(top_node),
            "is_solved": top_state.is_solved,
            "number_of_steps": top_state.max_transforms,
            "number_of_precursors": len(top_state.mols),
            "number_of_precursors_in_stock": sum(top_state.in_stock_list),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
            "policy_used_counts": dict(policy_used_counts),
        }

    def _all_nodes(self) -> Sequence[Node]:
        # This is to keep backwards compatibility, this should be investigate further
        if repr(self.scorer) == "state score":
            return list(self.search_tree.graph())
        return [node for node in self.search_tree.graph() if not node.children()]


class ReactionTree:
    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either FixedRetroReaction or UniqueMolecule objects.

    :ivar has_repeating_patterns: if the graph has repetitive elements
    :ivar graph: the bipartite graph
    :ivar is_solved: if all of the leaf nodes are in stock
    :ivar root: the root of the tree
    """

    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.root = none_molecule()
        self.has_repeating_patterns: bool = False
        self.is_solved: bool = False

    @classmethod
    def from_analysis(
        cls, analysis: TreeAnalysis, from_node: Node = None
    ) -> "ReactionTree":
        """
        Create a reaction from a tree analysis.

        The single route can be from a given node in the search tree if the `from_node`
        argument is given. If it is not given, the top scoring node is used.

        :param analysis: the analysis to base the reaction tree on
        :param from_node: the end node of the route, defaults to None
        :returns: the reaction tree
        """
        if not from_node:
            from_node = analysis.best_node()

        actions, nodes = analysis.search_tree.route_to_node(from_node)
        return ReactionTreeFromMcts(actions=actions, nodes=nodes).tree

    @classmethod
    def from_dict(cls, tree_dict: StrDict) -> "ReactionTree":
        """
        Create a new ReactionTree by parsing a dictionary.

        This is supposed to be the opposite of ``to_dict``,
        but because that format loses information, the returned
        object is not a full copy:
        * The stock will only contain the list of molecules marked as ``in_stock`` in the dictionary.
        * The reaction nodes will be of type `FixedRetroReaction`

        The returned object should be sufficient to e.g. generate an image of the route.

        :param tree_dict: the dictionary representation
        :returns: the reaction tree
        """
        return ReactionTreeFromDict(tree_dict).tree

    def depth(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> int:
        """
        Return the depth of a node in the route

        :param node: the query node
        :return: the depth
        """
        return self.graph.nodes[node].get("depth", -1)

    def distance_to(self, other: "ReactionTree", content: str = "both") -> float:
        """
        Calculate the distance to another reaction tree

        This is a tree edit distance, with unit cost to
        insert and deleted nodes, and the Jaccard distance for substituting nodes

        :param other: the reaction tree to compare to
        :param content: determine what part of the tree to include in the calculation
        :return: the distance between the routes
        """
        self_wrapper = ReactionTreeWrapper(self, content)
        other_wrapper = ReactionTreeWrapper(other, content)
        return self_wrapper.distance_to(other_wrapper)

    def in_stock(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> bool:
        """
        Return if a node in the route is in stock

        Note that is a property set on creation and as such is not updated.

        :param node: the query node
        :return: if the molecule is in stock
        """
        return self.graph.nodes[node].get("in_stock", False)

    def leafs(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecules nodes of the reaction tree that has no predecessors,
        i.e. molecules that has not been broken down

        :yield: the next leaf molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule) and not self.graph[node]:
                yield node

    def molecules(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecule nodes of the reaction tree

        :yield: the next molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule):
                yield node

    def reactions(self) -> Iterable[FixedRetroReaction]:
        """
        Generates the reaction nodes of the reaction tree

        :yield: the next reaction in the tree
        """
        for node in self.graph:
            if not isinstance(node, Molecule):
                yield node

    def to_dict(self) -> StrDict:
        """
        Returns the reaction tree as a dictionary in a pre-defined format.

        :return: the reaction tree
        """
        return self._build_dict(self.root)

    def to_image(
        self,
        in_stock_colors: FrameColors = None,
        show_all: bool = True,
    ) -> PilImage:
        """
        Return a pictoral representation of the route

        :raises ValueError: if image could not be produced
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :param show_all: if True, also show nodes that are marked as hidden
        :return: the image of the route
        """

        def show(node_):
            return not self.graph.nodes[node_].get("hide", False)

        in_stock_colors = in_stock_colors or {True: "green", False: "orange"}
        molecules = []
        frame_colors = []
        for node in self.molecules():
            if not show_all and not show(node):
                continue
            molecules.append(node)
            frame_colors.append(in_stock_colors[self.in_stock(node)])

        reactions = [node for node in self.reactions() if show_all or show(node)]
        edges = [
            (node1, node2)
            for node1, node2 in self.graph.edges()
            if show_all or (show(node1) and show(node2))
        ]

        try:
            return make_graphviz_image(molecules, reactions, edges, frame_colors)
        except FileNotFoundError as err:
            raise ValueError(str(err))

    def to_json(self) -> str:
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        """
        return json.dumps(self.to_dict(), sort_keys=False, indent=2)

    def _build_dict(
        self, node: Union[UniqueMolecule, FixedRetroReaction], dict_: StrDict = None
    ) -> StrDict:
        if dict_ is None:
            dict_ = {}

        dict_["type"] = "mol" if isinstance(node, Molecule) else "reaction"
        dict_["hide"] = self.graph.nodes[node].get("hide", False)
        dict_["smiles"] = node.smiles
        if isinstance(node, UniqueMolecule):
            dict_["is_chemical"] = True
            dict_["in_stock"] = self.in_stock(node)
        elif isinstance(node, FixedRetroReaction):
            dict_["is_reaction"] = True
            dict_["metadata"] = dict(node.metadata)
        else:
            raise ValueError(
                f"This is an invalid reaction tree. Unknown node type {type(node)}"
            )

        dict_["children"] = []

        for child in self.graph.successors(node):
            child_dict = self._build_dict(child)
            dict_["children"].append(child_dict)

        if not dict_["children"]:
            del dict_["children"]
        return dict_


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
        nodes, scores = analysis.sort_nodes(min_return=min_nodes)
        reaction_trees = [ReactionTree.from_analysis(analysis, node) for node in nodes]
        all_scores = [{repr(analysis.scorer): score} for score in scores]
        return cls(
            reaction_trees=reaction_trees,
            nodes=nodes,
            scores=scores,
            all_scores=all_scores,
        )

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
        self, n_clusters: int, max_clusters: int = 5, timeout: int = None, **kwargs: Any
    ) -> np.ndarray:
        """
        Cluster the route collection into a number of clusters.

        Additional arguments to the distance or clustering algorithm
        can be passed in as key-word arguments.

        If the number of reaction trees are less than 3, no clustering will be performed

        :param n_clusters: the desired number of clusters, if less than 2 triggers optimization
        :param max_clusters: the maximum number of clusters to consider
        :param timeout: if given, return no clusters if distance calculation is taking longer time
        :return: the cluster labels
        """
        if len(self.reaction_trees) < 3:
            return np.asarray([])
        content = kwargs.pop("content", "both")
        try:
            distances = self.distance_matrix(content=content, timeout=timeout)
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
        self,
        content: str = "both",
        recreate: bool = False,
        timeout: int = None,
    ) -> np.ndarray:
        """
        Compute the distance matrix between each pair of reaction trees

        :param content: determine what part of the tree to include in the calculation
        :param recreate: if False, use a cached one if available
        :param timeout: if given, raises an exception if timeout is taking longer time
        :return: the square distance matrix
        """
        if self._distance_matrix.get(content) is not None and not recreate:
            return self._distance_matrix[content]
        distances = np.zeros([len(self), len(self)])
        distance_wrappers = [
            ReactionTreeWrapper(rt, content) for rt in self.reaction_trees
        ]
        time0 = time.perf_counter()
        for i, iwrapper in enumerate(distance_wrappers):
            # fmt: off
            for j, jwrapper in enumerate(distance_wrappers[i + 1:], i + 1):
                distances[i, j] = iwrapper.distance_to(jwrapper)
                distances[j, i] = distances[i, j]
            # fmt: on
            time_past = time.perf_counter() - time0
            if timeout is not None and time_past > timeout:
                raise ValueError(f"Unable to compute distance matrix in {timeout} s")
        self._distance_matrix[content] = distances
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
