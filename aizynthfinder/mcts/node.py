""" Module containing a class that represents a node in the search tree.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np

from aizynthfinder.chem import RetroReaction, TreeMolecule
from aizynthfinder.mcts.state import State
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import StrDict, List, Optional
    from aizynthfinder.utils.serialization import (
        MoleculeDeserializer,
        MoleculeSerializer,
    )
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.mcts.mcts import SearchTree


class NodeUnexpectedBehaviourException(Exception):
    """Exception that is raised if the tree search is behaving unexpectedly."""


class Node:
    """
    A node in the search tree.

    The children are instantiated lazily for efficiency: only when
    a child is selected the reaction to create that child is applied.

    Properties of an instantiated children to a node can be access with:

    .. code-block::

        children_attr = node[child]

    the return value is a dictionary with keys "action", "value", "prior"
    and "visitations".

    :ivar is_expanded: if the node has had children added to it
    :ivar is_expandable: if the node is expandable
    :ivar parent: the parent node
    :ivar state: the internal state of the node
    :ivar tree: the tree owning this node

    :param state: the state of the node
    :param owner: the tree that owns this node
    :param config: settings of the tree search algorithm
    :param parent: the parent node, defaults to None
    """

    def __init__(
        self, state: State, owner: SearchTree, config: Configuration, parent=None
    ):
        self.state = state
        self._config = config
        self._expansion_policy = config.expansion_policy
        self._filter_policy = config.filter_policy
        self.tree = owner
        self.is_expanded: bool = False
        self.is_expandable: bool = not self.state.is_terminal
        self.parent = parent

        self._children_values: List[float] = []
        self._children_priors: List[float] = []
        self._children_visitations: List[int] = []
        self._children_actions: List[RetroReaction] = []
        self._children: List[Optional[Node]] = []

        self._logger = logger()

    def __getitem__(self, node: "Node") -> StrDict:
        idx = self._children.index(node)
        return {
            "action": self._children_actions[idx],
            "value": self._children_values[idx],
            "prior": self._children_priors[idx],
            "visitations": self._children_visitations[idx],
        }

    @classmethod
    def create_root(
        cls, smiles: str, tree: SearchTree, config: Configuration
    ) -> "Node":
        """
        Create a root node for a tree using a SMILES.

        :param smiles: the SMILES representation of the root state
        :param tree: the search tree
        :param config: settings of the tree search algorithm
        :return: the created node
        """
        mol = TreeMolecule(parent=None, transform=0, smiles=smiles)
        state = State(mols=[mol], config=config)
        return Node(state=state, owner=tree, config=config)

    @classmethod
    def from_dict(
        cls,
        dict_: StrDict,
        tree: SearchTree,
        config: Configuration,
        molecules: MoleculeDeserializer,
        parent: "Node" = None,
    ) -> "Node":
        """
        Create a new node from a dictionary, i.e. deserialization

        :param dict_: the serialized node
        :param tree: the search tree
        :param config: settings of the tree search algorithm
        :param molecules: the deserialized molecules
        :param parent: the parent node
        :return: a deserialized node
        """
        state = State.from_dict(dict_["state"], config, molecules)
        node = Node(state=state, owner=tree, config=config, parent=parent)
        node.is_expanded = dict_["is_expanded"]
        node.is_expandable = dict_["is_expandable"]
        node._children_values = dict_["children_values"]
        node._children_priors = dict_["children_priors"]
        node._children_visitations = dict_["children_visitations"]
        node._children_actions = []
        for action in dict_["children_actions"]:
            mol = molecules.get_tree_molecules([action["mol"]])[0]
            node._children_actions.append(
                RetroReaction(
                    mol,
                    action["smarts"],
                    action["index"],
                    action.get("metadata", {}),
                )
            )

        node._children = [
            Node.from_dict(child, tree, config, molecules, parent=node)
            if child
            else None
            for child in dict_["children"]
        ]
        return node

    def backpropagate(self, child: "Node", value_estimate: float) -> None:
        """
        Update the number of visitations of a particular child and its value.

        :param child: the child node
        :param value_estimate: the value to add to the child value
        """
        idx = self._children.index(child)
        self._children_visitations[idx] += 1
        self._children_values[idx] += value_estimate

    def children(self) -> List["Node"]:
        """
        Returns all of the instantiated children

        :return: the children
        """
        return [child for child in self._children if child]

    def children_view(self) -> StrDict:
        """
        Creates a view of the children attributes. Each of the
        list returned is a new list, although the actual children
        are not copied.

        The return dictionary will have keys "actions", "values",
        "priors", "visitations" and "objects".

        :return: the view
        """
        return {
            "actions": list(self._children_actions),
            "values": list(self._children_values),
            "priors": list(self._children_priors),
            "visitations": list(self._children_visitations),
            "objects": list(self._children),
        }

    def expand(self) -> None:
        """
        Expand the node.

        Expansion is the process of creating the children of the node,
        without instantiating a child object. The actions and priors are
        taken from the policy network.
        """
        if self.is_expanded:
            msg = f"Oh no! This node is already expanded. id={id(self)}"
            self._logger.debug(msg)
            raise NodeUnexpectedBehaviourException(msg)

        if self.is_expanded or not self.is_expandable:
            return

        self.is_expanded = True

        # Calculate the possible actions, fill the child_info lists
        # Actions by default only assumes 1 set of reactants
        (
            self._children_actions,
            self._children_priors,
        ) = self._expansion_policy.get_actions(self.state.mols)
        nactions = len(self._children_actions)
        self._children_visitations = [1] * nactions
        self._children = [None] * nactions
        if self._config.use_prior:
            self._children_values = list(self._children_priors)
        else:
            self._children_values = [self._config.default_prior] * nactions

    def is_terminal(self) -> bool:
        """
        Node is terminal if its unexpandable, or the internal state is terminal (solved)

        :return: the terminal attribute of the node
        """
        return not self.is_expandable or self.state.is_terminal

    def promising_child(self) -> Optional["Node"]:
        """
        Return the child with the currently highest Q+U

        The selected child will be instantiated if it has not been already.

        If no actions could be found that were applicable, the method will
        return None.

        :return: the child
        """

        scores = self._children_q() + self._children_u()
        indices = np.where(scores == scores.max())[0]
        index = np.random.choice(indices)

        child = self._select_child(index)
        if not child and max(self._children_values) > 0:
            return self.promising_child()

        if not child:
            self._logger.debug(
                "Returning None from promising_child() because there were no applicable action"
            )
            self.is_expanded = False
            self.is_expandable = False

        return child

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :return: the serialized node
        """
        return {
            "state": self.state.serialize(molecule_store),
            "children_values": [float(value) for value in self._children_values],
            "children_priors": [float(value) for value in self._children_priors],
            "children_visitations": self._children_visitations,
            "children_actions": [
                {
                    "mol": molecule_store[action.mol],
                    "smarts": action.smarts,
                    "index": action.index,
                    "metadata": dict(action.metadata),
                }
                for action in self._children_actions
            ],
            "children": [
                child.serialize(molecule_store) if child else None
                for child in self._children
            ],
            "is_expanded": self.is_expanded,
            "is_expandable": self.is_expandable,
        }

    def _check_child_reaction(self, reaction: RetroReaction) -> bool:
        if not reaction.reactants:
            self._logger.debug(
                f"Reaction {reaction.smarts} on {reaction.mol.smiles} did not produce any reactants"
            )
            return False

        # fmt: off
        reactants0 = reaction.reactants[0]
        if len(reaction.reactants) == 1 and len(reactants0) == 1 and reaction.mol == reactants0[0]:
            return False

        check_cycling = self._config.prune_cycles_in_search and reaction.mol.parent
        if check_cycling and len(reaction.reactants) == 1 and (reaction.mol.parent in reactants0):
            self._logger.debug(
                "Removing cycle-producing reaction on mol %s with parent %s and transformation %s"
                % (reaction.mol.smiles, reaction.mol.parent.smiles, reaction.smarts)
            )
            return False
        # fmt: on

        return True

    def _children_q(self) -> np.ndarray:
        return np.array(self._children_values) / np.array(self._children_visitations)

    def _children_u(self) -> np.ndarray:
        total_visits = np.log(np.sum(self._children_visitations))
        child_visits = np.array(self._children_visitations)
        return self._config.C * np.sqrt(2 * total_visits / child_visits)

    def _create_children_nodes(
        self, states: List[State], child_idx: int
    ) -> List[Optional["Node"]]:
        new_nodes = []
        first_child_idx = child_idx
        for state_index, state in enumerate(states):
            # If there's more than one outcome, the lists need be expanded
            if state_index > 0:
                child_idx = self._expand_children_lists(first_child_idx, state_index)

            if self._filter_child_reaction(self._children_actions[child_idx]):
                self._children_values[child_idx] = -1e6
            else:
                self._children[child_idx] = Node(
                    state=state, owner=self.tree, config=self._config, parent=self
                )
                new_nodes.append(self._children[child_idx])
        return new_nodes

    def _expand_children_lists(self, old_index: int, action_index: int) -> int:
        new_action = self._children_actions[old_index].copy(index=action_index)
        self._children_actions.append(new_action)
        self._children_priors.append(self._children_priors[old_index])
        self._children_values.append(self._children_values[old_index])
        self._children_visitations.append(self._children_visitations[old_index])
        self._children.append(None)
        return len(self._children) - 1

    def _filter_child_reaction(self, reaction: RetroReaction) -> bool:
        if not self._filter_policy.selection:
            return False
        feasible, prob = self._filter_policy.feasibility(reaction)
        if not feasible:
            self._logger.debug(
                f"Reaction {reaction.reaction_smiles()} from template {reaction.smarts } "
                f"was filtered out with prob {prob}"
            )
            return True

        self._logger.debug(
            f"Reaction {reaction.reaction_smiles()} from template {reaction.smarts } was kept with prob {prob}"
        )
        return False

    def _select_child(self, child_idx: int) -> Optional["Node"]:
        """
        Selecting a child node implies instantiating the children nodes

        The algorithm is:
        * If the child has already been instantiated, return immediately
        * Apply the reaction associated with the child
        * If the application of the action failed, set value to -1e6 and return None
        * Create a new state array, one new state for each of the reaction outcomes
        * Create new child nodes
            - If a filter policy is available and the reaction outcome is unlikely
              set value of child to -1e6
        * Select a random node of the feasible ones to return
        """
        if self._children[child_idx]:
            return self._children[child_idx]

        reaction = self._children_actions[child_idx]
        if not reaction.reactants:
            reaction.apply()

        if not self._check_child_reaction(reaction):
            self._children_values[child_idx] = -1e6
            return None

        keep_mols = [mol for mol in self.state.mols if mol is not reaction.mol]
        new_states = [
            State(keep_mols + list(reactants), self._config)
            for reactants in reaction.reactants
        ]
        new_nodes = self._create_children_nodes(new_states, child_idx)

        if new_nodes:
            return np.random.choice(new_nodes)
        return None
