""" Module containing a class that holds the tree search
"""
from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING

import networkx as nx

from aizynthfinder.chem import MoleculeDeserializer, MoleculeSerializer
from aizynthfinder.search.mcts.node import MctsNode, ParetoMctsNode
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import List, Optional, Sequence, Union


_MODE2NODECLASS = {
    "single-objective": MctsNode,
    "weighted-sum": MctsNode,
    "multi-objective": ParetoMctsNode,
}


class MctsSearchTree:
    """
    Encapsulation of the search tree.

    :ivar root: the root node
    :ivar config: the configuration of the search tree

    :param config: settings of the tree search algorithm
    :param root_smiles: the root will be set to a node representing this molecule, defaults to None
    """

    def __init__(
        self, config: Configuration, root_smiles: Optional[str] = None
    ) -> None:
        self._logger = logger()
        self.profiling = {
            "expansion_calls": 0,
            "reactants_generations": 0,
            "iterations": 0,
        }
        self.config = config
        self.mode = self._check_mode()
        self._logger.debug(f"MCTS mode: {self.mode}")

        if root_smiles:
            self.root: Optional[MctsNode] = _MODE2NODECLASS[self.mode].create_root(
                smiles=root_smiles, tree=self, config=config
            )
        else:
            self.root = None

        self._graph: Optional[nx.DiGraph] = None

        # For backward compatibility
        if "search_reward" in config.search.algorithm_config:
            config.search.algorithm_config["search_rewards"] = [
                config.search.algorithm_config["search_reward"]
            ]
        config_rewards = config.search.algorithm_config["search_rewards"]
        self._logger.setLevel(logging.WARNING)  # Supress logging from `make_subset`
        self._logger.debug(f"Selecting reward scorers: {config_rewards}")
        self.reward_scorer = self.config.scorers.make_subset(config_rewards)
        self._logger.setLevel(logging.DEBUG)
        if self.mode == "single-objective":
            self.reward_scorer_name = config_rewards[0]

    @classmethod
    def from_json(cls, filename: str, config: Configuration) -> "MctsSearchTree":
        """
        Create a new search tree by deserialization from a JSON file

        :param filename: the path to the JSON node
        :param config: the configuration of the search
        :return: a deserialized tree
        """
        tree = MctsSearchTree(config)
        with open(filename, "r") as fileobj:
            dict_ = json.load(fileobj)
        mol_deser = MoleculeDeserializer(dict_["molecules"])
        tree.root = _MODE2NODECLASS[tree.mode].from_dict(
            dict_["tree"], tree, config, mol_deser
        )
        return tree

    def backpropagate(self, from_node: MctsNode) -> None:
        """
        Backpropagate the value estimate and update all nodes from a
        given node all the way to the root.

        :param from_node: the end node of the route to update
        """
        value_estimate = self.compute_reward(from_node)

        current = from_node
        while current is not self.root:
            parent = current.parent
            # For mypy, parent should never by None unless current is the root
            assert parent is not None
            parent.backpropagate(current, value_estimate)  # type: ignore
            current = parent

    def compute_reward(self, node: MctsNode) -> Union[float, Sequence[float]]:
        """
        Compute the reward of a node in the search tree, using one of

            1. Single-objective
            2. Multi-objective
            3. Weighted-sum of multiple specified rewards

        :param node: the node to compute the reward for
        :returns: the value from the scorer(s)
        """
        if self.mode == "single-objective":
            return self.reward_scorer[self.reward_scorer_name](node)

        if self.mode == "multi-objective":
            return self.reward_scorer.score_vector(node)

        return self.reward_scorer.weighted_score(
            node, self.config.search.algorithm_config["search_rewards_weights"]
        )

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
            for grandchild in node.children:
                add_node(grandchild)

        self._graph = nx.DiGraph()
        # Always add the root
        self._graph.add_node(self.root)
        for child in self.root.children:
            add_node(child)
        return self._graph

    def nodes(self) -> List[MctsNode]:
        """Return all the nodes in the search tree"""
        return list(self.graph())

    def one_iteration(self) -> bool:
        """
        Perform one iteration of
            1. Selection
            2. Expansion
            3. Rollout
            4. Backpropagation

        :return: if a solution was found
        """
        self.profiling["iterations"] += 1
        leaf = self.select_leaf()
        leaf.expand()
        while not leaf.is_terminal():
            child = leaf.promising_child()
            if child:
                child.expand()
                leaf = child

        self.backpropagate(leaf)
        return leaf.state.is_solved

    def select_leaf(self) -> MctsNode:
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

    def _check_mode(self) -> str:
        # if no objective weights are supplied, use multi-objective search
        # if only one objective is specified, search will in be in multi objective mode,
        #  but it is simply a vectoried version of single objective mode
        nrewards = len(self.config.search.algorithm_config["search_rewards"])
        nweights = len(self.config.search.algorithm_config["search_rewards_weights"])

        if nrewards == 1:
            mode = "single-objective"
        elif nweights == 0:
            mode = "multi-objective"
        else:
            mode = "weighted-sum"

        if mode == "weighted-sum" and nrewards != nweights:
            raise ValueError(
                "number of reward weights are expected to be the same as number of objectives,"
                f"currently have {nweights} weights and {nrewards} objectives)"
            )
        return mode
