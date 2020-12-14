""" Module containing classes to perform analysis of the tree search results.
"""
import json
from collections import defaultdict

import numpy as np
import networkx as nx


from aizynthfinder.chem import (
    Molecule,
    hash_reactions,
)
from aizynthfinder.utils.image import make_graphviz_image
from aizynthfinder.utils.analysis_helpers import (
    ReactionTreeFromDict,
    ReactionTreeFromMcts,
    CombinedReactionTrees,
)
from aizynthfinder.utils.route_clustering import ReactionTreeWrapper, ClusteringHelper


class TreeAnalysis:
    """
    Class that encapsulate various analysis that can be
    performed on a search tree.

    :ivar scorer: the object used to score the nodes
    :type scorer: Scorer
    :ivar search_tree: the search tree
    :vartype search_tree: SearchTree

    :param scorer: the object used to score the nodes, defaults to StateScorer
    :type scorer: Scorer, optional
    :parameter search_tree: the search tree to do the analysis on
    :type search_tree: SearchTree
    """

    def __init__(self, search_tree, scorer=None):
        self.search_tree = search_tree
        if scorer is None:
            # Do import here to avoid circular imports
            from aizynthfinder.context.scoring import StateScorer

            self.scorer = StateScorer()
        else:
            self.scorer = scorer

    def best_node(self):
        """
        Returns the node with the highest score.
        If several nodes have the same score, it will return the first

        :return: the top scoring node
        :rtype: Node
        """
        nodes = self._all_nodes()
        sorted_nodes, _ = self.scorer.sort(nodes)
        return sorted_nodes[0]

    def sort_nodes(self, min_return=5, max_return=25):
        """
        Sort and select the nodes, so that the best scoring routes are returned.
        The algorithm filter away identical routes and returns at minimum the number specified.
        If multiple alternative routes have the same score as the n'th route, they will be included and returned.

        :param min_return: the minium number of routes to return, defaults to 5
        :type min_return: int, optional
        :param max_return: the maximum number of routes to return
        :type max_return: int, optional
        :return: the nodes
        :rtype: list of Node
        :return: the score
        :rtype: list of float
        """
        nodes = self._all_nodes()
        sorted_nodes, sorted_scores = self.scorer.sort(nodes)

        if len(nodes) <= min_return:
            return sorted_nodes, sorted_scores

        seen_hashes = set()
        best_nodes = []
        best_scores = []
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

    def tree_statistics(self):
        """
        Returns statiscs of the tree

        Currently it returns the number of nodes, the maximum number of transforms,
        maximum number of children, top score, if the top score route is solved,
        the number of molecule in the top score node, and information on pre-cursors

        :return: the statistics
        :rtype: dict
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

        policy_used_counts = defaultdict(lambda: 0)
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

    def _all_nodes(self):
        # This is to keep backwards compatibility, this should be investigate further
        if repr(self.scorer) == "state score":
            return list(self.search_tree.graph())
        return [node for node in self.search_tree.graph() if not node.children()]


class ReactionTree:
    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either FixedRetroReaction or UniqueMolecule objects.

    :ivar has_repeating_patterns: if the graph has repetetive elements
    :vartype has_repeating_patterns: bool
    :ivar graph: the bipartite graph
    :vartype graph: networkx.DiGraph
    :ivar is_solved: if all of the leaf nodes are in stock
    :vartype is_solved: bool
    :ivar root: the root of the tree
    :vartype root: UniqueMolecule
    """

    def __init__(self):
        self.graph = nx.DiGraph()
        self.root = None
        self.has_repeating_patterns = False
        self.is_solved = False

    @classmethod
    def from_analysis(cls, analysis, from_node=None):
        """
        Create a reaction from a tree analysis.

        The single route can be from a given node in the search tree if the `from_node`
        argument is given. If it is not given, the top scoring node is used.

        :param analysis: the analysis to base the reaction tree on
        :type analysis: TreeAnalysis
        :param from_node: the end node of the route, defaults to None
        :type from_node: Node, optional
        :returns: the reaction tree
        :rtype: ReactionTree
        """
        if not from_node:
            from_node = analysis.best_node()

        actions, nodes = analysis.search_tree.route_to_node(from_node)
        return ReactionTreeFromMcts(actions=actions, nodes=nodes).tree

    @classmethod
    def from_dict(cls, tree_dict):
        """
        Create a new ReactionTree by parsing a dictionary.

        This is supposed to be the opposite of ``to_dict``,
        but because that format loses information, the returned
        object is not a full copy:
        * The stock will only contain the list of molecules marked as ``in_stock`` in the dictionary.
        * The reaction nodes will be of type `FixedRetroReaction`

        The returned object should be sufficient to e.g. generate an image of the route.

        :param tree_dict: the dictionary representation
        :type tree_dict: dict
        :returns: the reaction tree
        :rtype: ReactionTree
        """
        return ReactionTreeFromDict(tree_dict).tree

    def depth(self, node):
        """
        Return the depth of a node in the route

        :param node: the query node
        :type node: UniqueMolecule or Reaction
        :return: the depth
        :rtype: float
        """
        return self.graph.nodes[node].get("depth", -1)

    def distance_to(self, other, content="both"):
        """
        Calculate the distance to another reaction tree

        This is a tree edit distance, with unit cost to
        insert and deleted nodes, and the Jaccard distance for substituting nodes

        :param other: the reaction tree to compare to
        :type other: ReactionTree
        :param content: determine what part of the tree to include in the calculation
        :type content: TreeContent or str, optional
        :return: the distance between the routes
        :rtype: float
        """
        self_wrapper = ReactionTreeWrapper(self, content)
        other_wrapper = ReactionTreeWrapper(other, content)
        return self_wrapper.distance_to(other_wrapper)

    def in_stock(self, node):
        """
        Return if a node in the route is in stock

        Note that is a property set on creation and as such is not updated.

        :param node: the query node
        :type node: UniqueMolecule
        :return: if the molecule is in stock
        :rtype: bool
        """
        return self.graph.nodes[node].get("in_stock", False)

    def leafs(self):
        """
        Generates the molecules nodes of the reaction tree that has no predecessors,
        i.e. molecules that has not been broken down

        :yield: the next leaf molecule in the tree
        :rtype: UniqueMolecule
        """
        for node in self.graph:
            if isinstance(node, Molecule) and not self.graph[node]:
                yield node

    def molecules(self):
        """
        Generates the molecule nodes of the reaction tree

        :yield: the next molecule in the tree
        :rtype: UniqueMolecule
        """
        for node in self.graph:
            if isinstance(node, Molecule):
                yield node

    def reactions(self):
        """
        Generates the reaction nodes of the reaction tree

        :yield: the next reaction in the tree
        :rtype: Reaction
        """
        for node in self.graph:
            if not isinstance(node, Molecule):
                yield node

    def to_dict(self):
        """
        Returns the reaction tree as a dictionary in a pre-defined format.

        :return: the reaction tree
        :rtype: dict
        """
        return self._build_dict(self.root)

    def to_image(self, in_stock_colors={True: "green", False: "orange"}, show_all=True):
        """
        Return a pictoral representation of the route

        :raises ValueError: if image could not be produced
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :type in_stock_colors: dict, optional
        :param show_all: if True, also show nodes that are marked as hidden
        :type show_all: bool, optional
        :return: the image of the route
        :rtype: PIL.Image
        """

        def show(node):
            return not self.graph.nodes[node].get("hide", False)

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

    def to_json(self):
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        :rtype: str
        """
        return json.dumps(self.to_dict(), sort_keys=False, indent=2)

    def _build_dict(self, node, dict_=None):
        if dict_ is None:
            dict_ = {}

        dict_["type"] = "mol" if isinstance(node, Molecule) else "reaction"
        dict_["hide"] = self.graph.nodes[node].get("hide", False)
        dict_["smiles"] = node.smiles
        if isinstance(node, Molecule):
            dict_["is_chemical"] = True
            dict_["in_stock"] = self.in_stock(node)
        else:
            dict_["is_reaction"] = True
            dict_["metadata"] = dict(node.metadata)

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
    :vartype all_scores: list of dict
    :ivar nodes: the top-ranked nodes
    :vartype nodes: list of Node
    :ivar scores: initial scores of top-ranked nodes
    :vartype scores: list of float
    :ivar reaction_trees: the reaction trees created from the top-ranked nodes
    :vartype reaction_trees: list of Reaction Tree
    :ivar clusters: the created clusters from the collection
    :vartype clusters: list of RouteCollection

    :param reaction_trees: the trees to base the collection on
    :type reaction_trees: list of ReactionTree
    """

    def __init__(self, reaction_trees, **kwargs):
        self._routes = [{} for _ in range(len(reaction_trees))]
        self.reaction_trees = reaction_trees
        self._update_route_dict(reaction_trees, "reaction_tree")

        self.nodes = self._unpack_kwarg("nodes", **kwargs)
        if not self.nodes:
            self.nodes = [None] * len(self.reaction_trees)

        self.scores = self._unpack_kwarg("scores", **kwargs)
        if not self.scores:
            self.scores = [np.nan] * len(self.reaction_trees)

        self.all_scores = self._unpack_kwarg("all_scores", **kwargs)
        if not self.all_scores:
            self.all_scores = [dict() for _ in range(len(self.reaction_trees))]

        self._dicts = self._unpack_kwarg("dicts", **kwargs)
        self._images = self._unpack_kwarg("images", **kwargs)
        self._jsons = self._unpack_kwarg("jsons", **kwargs)
        self._clusters = self._unpack_kwarg("clusters", **kwargs)
        self._distance_matrix = {}
        self._combined_reaction_trees = None

    @classmethod
    def from_analysis(cls, analysis, min_nodes):
        """
        Create a collection from a tree analysis.

        :param analysis: the tree analysis to use
        :type analysis: TreeAnalysis
        :param min_nodes: the minimum number of top-ranked nodes to consider
        :type min_nodes: int
        :return: the created collection
        :rtype: RouteCollection
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

    def __getitem__(self, index):
        if index < 0 or index >= len(self):
            raise IndexError("Index out of range")
        return self._routes[index]

    def __len__(self):
        return len(self.nodes)

    @property
    def dicts(self):
        """ Returns a list of dictionary representation of the routes
        """
        if self._dicts is None:
            self.make_dicts()
        return self._dicts

    @property
    def images(self):
        """ Returns a list of pictoral representation of the routes
        """
        if self._images is None:
            self.make_images()
        return self._images

    @property
    def jsons(self):
        """ Returns a list of JSON string representation of the routes
        """
        if self._jsons is None:
            self.make_jsons()
        return self._jsons

    def cluster(self, n_clusters, max_clusters=5, **kwargs):
        """
        Cluster the route collection into a number of clusters.

        Additional arguments to the distance or clustering algorithm
        can be passed in as key-word arguments.

        If the number of reaction trees are less than 3, no clustering will be performed

        :param n_clusters: the desired number of clusters, if less than 2 triggers optimization
        :type n_clusters: int
        :param max_clusters: the maximum number of clusters to consider
        :type max_clusters: int, optional
        :return: the cluster labels
        :rtype: np.ndarray
        """
        if len(self.reaction_trees) < 3:
            return []
        content = kwargs.pop("content", "both")
        labels = ClusteringHelper.cluster(
            self.distance_matrix(content=content),
            n_clusters,
            max_clusters=max_clusters,
            **kwargs
        )
        self._make_clusters(labels)
        return labels

    def combined_reaction_trees(self, recreate=False):
        """
        Return an object that combines all the reaction tree into a single reaction tree graph

        :param recreate: if False will return a cached object if available, defaults to False
        :type recreate: bool, optional
        :return: the combined trees
        :rtype: CombinedReactionTrees
        """
        if not self._combined_reaction_trees or recreate:
            self._combined_reaction_trees = CombinedReactionTrees(self.reaction_trees)
        return self._combined_reaction_trees

    def compute_scores(self, *scorers):
        """
        Compute new scores for all routes in this collection.
        They can then be accessed with the ``all_scores`` attribute.
        """
        if self.nodes[0]:
            list_ = self.nodes
        else:
            list_ = self.reaction_trees

        for scorer in scorers:
            for idx, score in enumerate(scorer(list_)):
                self.all_scores[idx][repr(scorer)] = score
        self._update_route_dict(self.all_scores, "all_score")

    def dict_with_scores(self):
        """
        Return the routes as dictionaries with all scores added
        to the root (target) node.

        :return: the routes as dictionaries
        :rtype: list of dict
        """
        dicts = []
        for dict_, scores in zip(self.dicts, self.all_scores):
            dicts.append(dict(dict_))
            dicts[-1]["scores"] = dict(scores)
        return dicts

    def distance_matrix(self, content="both", recreate=False):
        """
        Compute the distance matrix between each pair of reaction trees

        :param content: determine what part of the tree to include in the calculation
        :type content: TreeContent or str, optional
        :param recreate: if False, use a cached one if available
        :type recreate: bool, optional
        :return: the square distance matrix
        :rtype: numpy.ndarray
        """
        if self._distance_matrix.get(content) is not None and not recreate:
            return self._distance_matrix[content]
        distances = np.zeros([len(self), len(self)])
        distance_wrappers = [
            ReactionTreeWrapper(rt, content) for rt in self.reaction_trees
        ]
        for i, iwrapper in enumerate(distance_wrappers):
            # fmt: off
            for j, jwrapper in enumerate(distance_wrappers[i + 1:], i + 1):
                distances[i, j] = iwrapper.distance_to(jwrapper)
                distances[j, i] = distances[i, j]
            # fmt: on
        self._distance_matrix[content] = distances
        return distances

    def make_dicts(self):
        """ Convert all reaction trees to dictionaries
        """
        self._dicts = [tree.to_dict() for tree in self.reaction_trees]
        self._update_route_dict(self._dicts, "dict")

    def make_images(self):
        """ Convert all reaction trees to images
        """

        self._images = []
        for tree in self.reaction_trees:
            try:
                img = tree.to_image()
            except ValueError:
                self._images.append(None)
            else:
                self._images.append(img)
        self._update_route_dict(self._images, "image")

    def make_jsons(self):
        """ Convert all reaction trees to JSON strings
        """
        self._jsons = [tree.to_json() for tree in self.reaction_trees]
        self._update_route_dict(self._jsons, "json")

    def rescore(self, scorer):
        """
        Rescore the routes in the collection, and thereby re-order them.

        This will replace the ``scores`` attribute, and update the ``all_scores``
        attribute with another entry.

        :param scorer: the scorer to use
        :type scorer: Scorer
        """
        if self.nodes[0]:
            self.nodes, self.scores, sortidx = scorer.sort(
                self.nodes, return_sort_indices=True
            )
            self.reaction_trees = [self.reaction_trees[idx] for idx in sortidx]
        else:
            self.reaction_trees, self.scores, sortidx = scorer.sort(
                self.reaction_trees, return_sort_indices=True
            )
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

    def _make_clusters(self, clusters):
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

    def _select_subset(self, arr, selection):
        return [item for sel, item in zip(selection, arr) if sel]

    def _unpack_kwarg(self, key, **kwargs):
        if key not in kwargs:
            return None
        arr = kwargs[key]
        self._update_route_dict(arr, key[:-1])
        return arr

    def _update_route_dict(self, arr, key):
        for i, value in enumerate(arr):
            self._routes[i][key] = value
