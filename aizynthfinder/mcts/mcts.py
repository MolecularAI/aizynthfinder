""" Module containing a class that holds the tree search
"""
import json
import networkx as nx

from aizynthfinder.mcts.node import Node
from aizynthfinder.utils.serialization import MoleculeSerializer, MoleculeDeserializer


class SearchTree:
    """
    Encapsulation of the search tree.

    :ivar root: the root node
    :vartype root: Node

    :param config: settings of the tree search algorithm
    :type config: Configuration
    :param root_smiles: the root will be set to a node representing this molecule, defaults to None
    :type root_smiles: string, optional
    """

    def __init__(self, config, root_smiles=None):
        if root_smiles:
            self.root = Node.create_root(smiles=root_smiles, tree=self, config=config)
        else:
            self.root = None
        self._config = config
        self._graph = None

    @classmethod
    def from_json(cls, filename, config):
        """
        Create a new search tree by deserialization from a JSON file

        :param filename: the path to the JSON node
        :type dict_: str
        :return: a deserialized tree
        :rtype: SearchTree
        """
        tree = SearchTree(config)
        with open(filename, "r") as fileobj:
            dict_ = json.load(fileobj)
        mol_deser = MoleculeDeserializer(dict_["molecules"])
        tree.root = Node.from_dict(dict_["tree"], tree, config, mol_deser)
        return tree

    def backpropagate(self, from_node, value_estimate):
        """
        Backpropagate the value estimate and update all nodes from a
        given node all the way to the root.

        :param from_node: the end node of the route to update
        :type from_node: Node
        :param value_estimate: The score value to backpropagate
        :type value_estimate: float
        """
        current = from_node
        while current is not self.root:
            parent = current.parent
            parent.backpropagate(current, value_estimate)
            current = parent

    def graph(self, recreate=False):
        """
        Construct a directed graph object with the nodes as
        vertices and the actions as edges attribute "action".

        :param recreate: if True will construct the graph even though it is cached, defaults to False
        :type recreate: bool, optional
        :return: the graph object
        :rtype: networkx.DiGraph
        """
        if not recreate and self._graph:
            return self._graph

        def add_node(node):
            self._graph.add_edge(node.parent, node, action=node.parent[node]["action"])
            for child in node.children():
                add_node(child)

        self._graph = nx.DiGraph()
        # Always add the root
        self._graph.add_node(self.root)
        for child in self.root.children():
            add_node(child)
        return self._graph

    def select_leaf(self):
        """
        Traverse the tree selecting the most promising child at
        each step until leaf node returned.

        :return: the leaf node
        :rtype: Node
        """

        current = self.root
        while current.is_expanded and not current.state.is_solved:
            promising_child = current.promising_child()
            # If promising_child returns None it means that the node
            # is unexpandable, and hence we should break the loop
            if promising_child:
                current = promising_child
        return current

    def serialize(self, filename):
        """
        Seralize the search tree to a JSON file

        :param filename: the path to the JSON file
        :type filename: str
        """
        mol_ser = MoleculeSerializer()
        dict_ = {"tree": self.root.serialize(mol_ser), "molecules": mol_ser.store}
        with open(filename, "w") as fileobj:
            json.dump(dict_, fileobj, indent=2)

    @staticmethod
    def route_to_node(from_node):
        """
        Return the route to a give node to the root.

        Will return both the actions taken to go between the nodes,
        and the nodes in the route themselves.

        :param from_node: the end of the route
        :type from_node: Node
        :return: the route
        :rtype: tuple (list of action, list of Node)
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
