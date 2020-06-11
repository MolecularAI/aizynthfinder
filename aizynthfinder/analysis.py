""" Module containing classes to perform analysis of the tree search results.
"""
import json
import hashlib

import numpy as np
import networkx as nx

from aizynthfinder.chem import Molecule, UniqueMolecule
from aizynthfinder.utils.image import GraphvizReactionGraph


class TreeAnalysis:
    """
    Class that encapsulate various analysis that can be
    performed on a search tree.

    :ivar search_tree: the search tree
    :vartype search_tree: SearchTree

    :parameter search_tree: the search tree to do the analysis on
    :type search_tree: SearchTree
    """

    def __init__(self, search_tree):
        self.search_tree = search_tree

    def best_node(self):
        """
        Returns the node with the highest score.
        If several nodes have the same score, it will return the first

        :return: the top scoring node
        :rtype: Node
        """
        nodes = list(self.search_tree.graph())
        scores = [node.state.score for node in nodes]
        idx = np.argmax(scores)
        return nodes[idx]

    def sort_nodes(self, min_return=5):
        """
        Sort and select the nodes, so that the best scoring routes are returned.
        The algorithm filter away identical routes and returns at minimum the number specified.
        If multiple alternative routes have the same score as the n'th route, they will be included and returned.

        :param min_return: the minium number of routes to return, defaults to 5
        :type min_return: int, optional
        :return: the nodes
        :rtype: list of Node
        """
        nodes = list(self.search_tree.graph())
        scores = np.array([node.state.score for node in nodes])

        sortidx = np.argsort(scores)
        sorted_scores = scores[sortidx][::-1]
        sorted_nodes = [nodes[idx] for idx in sortidx][::-1]

        seen_hashes = set()
        best_nodes = []
        last_score = 1e16
        for score, node in zip(sorted_scores, sorted_nodes):
            if len(best_nodes) >= min_return and score < last_score:
                break
            route_actions, _ = self.search_tree.route_to_node(node)
            route_hash = self._routehash(route_actions)

            if route_hash in seen_hashes:
                continue
            seen_hashes.add(route_hash)
            best_nodes.append(node)
            last_score = score

        return best_nodes

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
        return {
            "number_of_nodes": len(nodes),
            "max_transforms": max(node.state.max_transforms for node in nodes),
            "max_children": max(len(node.children()) for node in nodes),
            "top_score": top_state.score,
            "is_solved": top_state.is_solved,
            "number_of_steps": top_state.max_transforms,
            "number_of_precursors": len(top_state.mols),
            "number_of_precursors_in_stock": sum(top_state.in_stock_list),
            "precursors_in_stock": mols_in_stock,
            "precursors_not_in_stock": mols_not_in_stock,
        }

    @staticmethod
    def _routehash(actions):
        reaction_str = ">>".join(action.reaction_smiles() for action in actions)
        reaction_list = reaction_str.replace(".", ">>").split(">>")
        hash_list = [
            hashlib.sha224(reactant.encode("utf8")).hexdigest()
            for reactant in reaction_list
        ]
        hash_list.sort()
        hash_list = ".".join(hash_list)
        return hashlib.sha224(hash_list.encode("utf8")).hexdigest()


class ReactionTree:
    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either Reaction objects and UniqueMolecule objects.

    :ivar graph: the bipartite graph
    :vartype graph: networkx.DiGraph
    :ivar root: the root of the tree
    :vartype root: UniqueMolecule
    """

    def __init__(self):
        self._stock = None
        self.graph = nx.DiGraph()
        self.root = None

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
        obj = cls()
        if not from_node:
            from_node = analysis.best_node()

        obj._stock = from_node.state.stock
        actions, nodes = analysis.search_tree.route_to_node(from_node)
        unique_mols = {}

        root_mol = nodes[0].state.mols[0]
        obj.root = root_mol.make_unique()
        unique_mols[id(root_mol)] = obj.root
        obj.graph.add_node(obj.root)

        for child, action in zip(nodes[1:], actions):
            obj._add_bipartite(child, action, unique_mols)

        return obj

    @classmethod
    def from_dict(cls, tree_dict):
        """
        Create a new ReactionTree by parsing a dictionary.

        This is supposed to be the opposite of ``to_dict``,
        but because that format loses information, the returned
        object is not a full copy:
        * The stock will only contain the list of molecules marked as ``in_stock`` in the dictionary.
        * The reaction nodes will not be of type Reaction, but a proxy class that only stores the SMILES

        The returned object should be sufficient to e.g. generate an image of the route.

        :param tree_dict: the dictionary representation
        :type tree_dict: dict
        :returns: the reaction tree
        :rtype: ReactionTree
        """

        class ProxyReaction:
            def __init__(self, smiles, metadata):
                self.smiles = smiles
                self.metadata = metadata

        obj = cls()
        obj._stock = []
        ReactionTree._parse_tree_dict(tree_dict, obj, ProxyReaction)

        return obj

    def to_dict(self):
        """
        Returns the reaction tree as a dictionary in a pre-defined format.

        :return: the reaction tree
        :rtype: dict
        """
        return self._build_dict(self.root)

    def to_image(self, in_stock_colors={True: "green", False: "orange"}):
        """
        Return a pictoral representation of the route

        :raises ValueError: if image could not be produced
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :type in_stock_colors: dict, optional
        :return: the image of the route
        :rtype: PIL.Image
        """
        img_graph = GraphvizReactionGraph()

        for node in self.graph:
            if isinstance(node, Molecule):
                img_graph.add_molecule(
                    node, frame_color=in_stock_colors[node in self._stock]
                )
            else:
                img_graph.add_reaction(node)
        for node1, node2 in self.graph.edges:
            img_graph.add_edge(node1, node2)

        try:
            return img_graph.to_image()
        except FileNotFoundError as err:
            raise ValueError(str(err))

    def to_json(self):
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        :rtype: str
        """
        return json.dumps(self.to_dict(), sort_keys=False, indent=2)

    def _add_bipartite(self, child, action, unique_mols):
        def unique_mol(molecule):
            id_ = id(molecule)
            if id_ not in unique_mols:
                unique_mols[id_] = molecule.make_unique()
            return unique_mols[id_]

        self.graph.add_edge(unique_mol(action.mol), action)
        for mol in child.state.mols:
            if mol.parent is action.mol:
                self.graph.add_edge(action, unique_mol(mol))

    def _build_dict(self, node, dict_=None):
        if dict_ is None:
            dict_ = {}

        dict_["type"] = "mol" if isinstance(node, Molecule) else "reaction"
        if isinstance(node, Molecule):
            dict_["smiles"] = node.smiles
            dict_["is_chemical"] = True
            dict_["in_stock"] = node in self._stock
        else:
            dict_["smiles"] = node.smiles
            dict_["is_reaction"] = True
            dict_["metadata"] = dict(node.metadata)

        dict_["children"] = []

        for child in self.graph.successors(node):
            child_dict = self._build_dict(child)
            dict_["children"].append(child_dict)

        if not dict_["children"]:
            del dict_["children"]
        return dict_

    @staticmethod
    def _parse_tree_dict(tree_dict, rt_object, reaction_cls):
        if tree_dict["type"] == "mol":
            node = UniqueMolecule(smiles=tree_dict["smiles"])
            if not rt_object.root:
                rt_object.root = node
            if tree_dict["in_stock"]:
                rt_object._stock.append(node)
        else:
            node = reaction_cls(
                smiles=tree_dict["smiles"], metadata=tree_dict.get("metadata", {})
            )
        rt_object.graph.add_node(node)

        for child_dict in tree_dict.get("children", []):
            child = ReactionTree._parse_tree_dict(child_dict, rt_object, reaction_cls)
            rt_object.graph.add_edge(node, child)
        return node


class RouteCollection:
    """
    Holds a collections of reaction routes, i.e. the top
    scored nodes, their scores and the reaction trees created
    from them.

    The class has the functionality to compute collective results
    for the different routes such as images.

    Properties of individual route can be obtained with simple indexing.

    .. code-block::

        route0 = collection[0]

    :ivar nodes: the top-ranked nodes
    :vartype nodes: list of Node
    :ivar scores: scores of top-ranked nodes
    :vartype scores: list of float
    :ivar reaction_trees: the reaction trees created from the top-ranked nodes
    :vartype reaction_trees: list of Reaction Tree

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
            self.score = [np.nan] * len(self.reaction_trees)

        self._dicts = self._unpack_kwarg("dicts", **kwargs)
        self._images = self._unpack_kwarg("images", **kwargs)
        self._jsons = self._unpack_kwarg("jsons", **kwargs)

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
        nodes = analysis.sort_nodes(min_return=min_nodes)
        scores = [node.state.score for node in nodes]
        reaction_trees = [ReactionTree.from_analysis(analysis, node) for node in nodes]
        return cls(reaction_trees=reaction_trees, nodes=nodes, scores=scores)

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

    def make_dicts(self):
        """ Convert all reaction trees to dictionaries
        """
        self._dicts = [tree.to_dict() for tree in self.reaction_trees]
        self._update_route_dict(self._dicts, "dict")

    def make_images(self):
        """ Convert all reaction trees to images
        """
        self._images = [tree.to_image() for tree in self.reaction_trees]
        self._update_route_dict(self._images, "image")

    def make_jsons(self):
        """ Convert all reaction trees to JSON strings
        """
        self._jsons = [tree.to_json() for tree in self.reaction_trees]
        self._update_route_dict(self._jsons, "json")

    def _unpack_kwarg(self, key, **kwargs):
        if key not in kwargs:
            return None
        arr = kwargs[key]
        self._update_route_dict(arr, key[:-1])
        return arr

    def _update_route_dict(self, arr, key):
        for i, value in enumerate(arr):
            self._routes[i][key] = value
