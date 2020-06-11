""" Module containing a class that represents a node in the search tree.
"""
import numpy as np

from aizynthfinder.chem import Reaction, TreeMolecule
from aizynthfinder.mcts.state import State
from aizynthfinder.utils.logging import logger


class NodeUnexpectedBehaviourException(Exception):
    """ Exception that is raised if the tree search is behaving unexpectedly.
    """


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
    :vartype is_expanded: bool
    :ivar parent: the parent node
    :ivar is_expandable: if the node is expandable
    :vartype is_expandable: bool
    :vartype parent: Node
    :ivar state: the internal state of the node
    :vartype state: State

    :param state: the state of the node
    :type state: State
    :param owner: the tree that owns this node
    :type owner: SearchTree
    :param config: settings of the tree search algorithm
    :type config: Configuration
    :param parent: the parent node, defaults to None
    :type parent: Node, optional
    """

    def __init__(self, state, owner, config, parent=None):
        self.state = state
        self._config = config
        self._policy = config.policy
        self._tree = owner
        self.is_expanded = False
        self.is_expandable = not self.state.is_terminal
        self.parent = parent

        self._children_values = []
        self._children_priors = []
        self._children_visitations = []
        self._children_actions = []
        self._children = []

        self._logger = logger()

    def __getitem__(self, node):
        idx = self._children.index(node)
        return {
            "action": self._children_actions[idx],
            "value": self._children_values[idx],
            "prior": self._children_priors[idx],
            "visitations": self._children_visitations[idx],
        }

    @classmethod
    def create_root(cls, smiles, tree, config):
        """
        Create a root node for a tree using a SMILES.

        :param smiles: the SMILES representation of the root state
        :type smiles: str
        :param tree: the search tree
        :type tree: SearchTree
        :param config: settings of the tree search algorithm
        :type config: Configuration
        :return: the created node
        :rtype: Node
        """
        mol = TreeMolecule(parent=None, transform=0, smiles=smiles)
        state = State(mols=[mol], config=config)
        return Node(state=state, owner=tree, config=config)

    @classmethod
    def from_dict(cls, dict_, tree, config, molecules, parent=None):
        """
        Create a new node from a dictionary, i.e. deserialization

        :param dict_: the serialized node
        :type dict_: dict
        :param tree: the search tree
        :type tree: SearchTree
        :param config: settings of the tree search algorithm
        :type config: Configuration
        :param molecules: the deserialized molecules
        :type molecules: MoleculeDeserializer
        :param parent: the parent node
        :type parent: Node
        :return: a deserialized node
        :rtype: Node
        """
        state = State.from_dict(dict_["state"], config, molecules)
        node = Node(state=state, owner=tree, config=config, parent=parent)
        node.is_expanded = dict_["is_expanded"]
        node.is_expandable = dict_["is_expandable"]
        node._children_values = dict_["children_values"]
        node._children_priors = dict_["children_priors"]
        node._children_visitations = dict_["children_visitations"]
        node._children_actions = [
            Reaction(molecules[action["mol"]], action["smarts"], action["index"])
            for action in dict_["children_actions"]
        ]
        node._children = [
            Node.from_dict(child, tree, config, molecules, parent=node)
            if child
            else None
            for child in dict_["children"]
        ]
        return node

    def backpropagate(self, child, value_estimate):
        """
        Update the number of visitations of a particular child and its value.

        :param child: the child node
        :type child: Node
        :param value_estimate: the value to add to the child value
        :type value_estimate: float
        """
        idx = self._children.index(child)
        self._children_visitations[idx] += 1
        self._children_values[idx] += value_estimate

    def children(self):
        """
        Returns all of the instantiated children

        :return: the children
        :rtype: list of Node
        """
        return [child for child in self._children if child]

    def children_view(self):
        """
        Creates a view of the children attributes. Each of the
        list returned is a new list, although the actual children
        are not copied.

        The return dictionary will have keys "actions", "values",
        "priors", "visitations" and "objects".

        :return: the view
        :rtype: dictionary of list
        """
        return {
            "actions": list(self._children_actions),
            "values": list(self._children_values),
            "priors": list(self._children_priors),
            "visitations": list(self._children_visitations),
            "objects": list(self._children),
        }

    def expand(self):
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
        self._children_actions, self._children_priors = self._policy.get_actions(
            self.state.mols
        )
        nactions = len(self._children_actions)
        self._children_visitations = [1] * nactions
        self._children = [None] * nactions
        if self._config.use_prior:
            self._children_values = list(self._children_priors)
        else:
            self._children_values = [self._config.default_prior] * nactions

    def is_terminal(self):
        """
        Node is terminal if its unexpandable, or the internal state is terminal (solved)

        :return: the terminal attribute of the node
        :rtype: bool
        """
        return not self.is_expandable or self.state.is_terminal

    def promising_child(self):
        """
        Return the child with the currently highest Q+U

        The selected child will be instantiated if it has not been already.

        If no actions could be found that were applicable, the method will
        return None.

        :return: the child
        :rtype: Node
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

    def serialize(self, molecule_store):
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :type molecule_store: MolecularSerializer
        :return: the serialized node
        :rtype: dict
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

    def _children_q(self):
        return np.array(self._children_values) / np.array(self._children_visitations)

    def _children_u(self):
        total_visits = np.log(np.sum(self._children_visitations))
        child_visits = np.array(self._children_visitations)
        return self._config.C * np.sqrt(2 * total_visits / child_visits)

    def _make_child_states(self, reaction):
        mols = [mol for mol in self.state.mols if mol is not reaction.mol]
        reactants_list = reaction.apply()

        if not reactants_list:
            self._logger.debug(
                "Reactants_list empty %s, for mol %s and transformation %s"
                % (repr(reactants_list), reaction.mol.smiles, reaction.smarts)
            )
            return []

        # fmt: off
        if (len(reactants_list) == 1 and len(reactants_list[0]) == 1 and reaction.mol == reactants_list[0][0]):
            return []
        # fmt: on

        return [
            State(mols + list(reactants), self._config) for reactants in reactants_list
        ]

    def _select_child(self, idx):
        reaction = self._children_actions[idx]

        if self._children[idx]:
            return self._children[idx]

        states = self._make_child_states(reaction)
        if not states:
            self._children_values[idx] = -1e6
            return None

        children = []
        for i, state in enumerate(states):
            child = Node(
                state=state, owner=self._tree, config=self._config, parent=self
            )
            children.append(child)
            # If there's more than one outcome, the lists need be expanded
            if i > 0:
                new_action = Reaction(
                    reaction.mol, reaction.smarts, index=i, metadata=reaction.metadata
                )
                self._children_actions.append(new_action)
                self._children_priors.append(self._children_priors[idx])
                self._children_values.append(self._children_values[idx])
                self._children_visitations.append(self._children_visitations[idx])
                self._children.append(child)
            else:
                self._children[idx] = child

        return np.random.choice(children)
