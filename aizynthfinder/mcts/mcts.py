""" Module containing a class that holds the tree search
"""
from __future__ import annotations
import json
from typing import TYPE_CHECKING

import networkx as nx

from aizynthfinder.mcts.node import Node
from aizynthfinder.utils.serialization import MoleculeSerializer, MoleculeDeserializer

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Tuple, List, Optional
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem import RetroReaction


class SearchTree:
    """
    Encapsulation of the search tree.

    :ivar root: the root node
    :ivar config: the configuration of the search tree

    :param config: settings of the tree search algorithm
    :param root_smiles: the root will be set to a node representing this molecule, defaults to None
    """

    def __init__(self, config: Configuration, root_smiles: str = None) -> None:
        if root_smiles:
            self.root: Optional[Node] = Node.create_root(
                smiles=root_smiles, tree=self, config=config
            )
        else:
            self.root = None
        self.config = config
        self._graph: Optional[nx.DiGraph] = None

    @classmethod
    def from_json(cls, filename: str, config: Configuration) -> "SearchTree":
        """
        Create a new search tree by deserialization from a JSON file

        :param filename: the path to the JSON node
        :param config: the configuration of the search
        :return: a deserialized tree
        """
        tree = SearchTree(config)
        with open(filename, "r") as fileobj:
            dict_ = json.load(fileobj)
        mol_deser = MoleculeDeserializer(dict_["molecules"])
        tree.root = Node.from_dict(dict_["tree"], tree, config, mol_deser)
        return tree

    def backpropagate(self, from_node: Node, value_estimate: float) -> None:
        """
        Backpropagate the value estimate and update all nodes from a
        given node all the way to the root.

        :param from_node: the end node of the route to update
        :param value_estimate: The score value to backpropagate
        """
        current = from_node
        while current is not self.root:
            parent = current.parent
            parent.backpropagate(current, value_estimate)
            current = parent

    def graph(self, recreate: bool = False) -> nx.DiGraph:
        """
        Construct a directed graph object with the nodes as
        vertices and the actions as edges attribute "action".

        :param recreate: if True will construct the graph even though it is cached, defaults to False
        :return: the graph object
        :raises ValueError: if the tree is not defined
        """
        if not self.root:
            raise ValueError("Root of search tree is not defined ")

        if not recreate and self._graph:
            return self._graph

        def add_node(node):
            self._graph.add_edge(node.parent, node, action=node.parent[node]["action"])
            for grandchild in node.children():
                add_node(grandchild)

        self._graph = nx.DiGraph()
        # Always add the root
        self._graph.add_node(self.root)
        for child in self.root.children():
            add_node(child)
        return self._graph

    def one_iteration(self) -> bool:
        """
        Perform one iteration of
            1. Selection
            2. Expansion
            3. Rollout
            4. Backpropagation

        :return: if a solution was found
        """
        leaf = self.select_leaf()
        leaf.expand()
        rollout_child = None
        while not leaf.is_terminal():
            child = leaf.promising_child()
            if not rollout_child:
                rollout_child = child
            if child:
                child.expand()
                leaf = child
        self.backpropagate(leaf, leaf.state.score)
        return leaf.state.is_solved

    def select_leaf(self) -> Node:
        """
        Traverse the tree selecting the most promising child at
        each step until leaf node returned.

        :return: the leaf node
        :raises ValueError: if the tree is not defined
        """
        if not self.root:
            raise ValueError("Root of search tree is not defined ")

        current = self.root
        while current.is_expanded and not current.state.is_solved:
            promising_child = current.promising_child()
            # If promising_child returns None it means that the node
            # is unexpandable, and hence we should break the loop
            if promising_child:
                current = promising_child
        return current

    def serialize(self, filename: str) -> None:
        """
        Serialize the search tree to a JSON file

        :param filename: the path to the JSON file
        :raises ValueError: if the tree is not defined
        """
        if not self.root:
            raise ValueError("Root of search tree is not defined ")

        mol_ser = MoleculeSerializer()
        dict_ = {"tree": self.root.serialize(mol_ser), "molecules": mol_ser.store}
        with open(filename, "w") as fileobj:
            json.dump(dict_, fileobj, indent=2)

    @staticmethod
    def route_to_node(from_node: Node) -> Tuple[List[RetroReaction], List[Node]]:
        """
        Return the route to a give node to the root.

        Will return both the actions taken to go between the nodes,
        and the nodes in the route themselves.

        :param from_node: the end of the route
        :return: the route
        """
        actions = []
        nodes = []
        current = from_node

        while current is not None:
            parent = current.parent
            if parent is not None:
                action = parent[current]["action"]
                actions.append(action)
            nodes.append(current)
            current = parent
        return actions[::-1], nodes[::-1]
