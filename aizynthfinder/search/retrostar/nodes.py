""" Module containing a classes representation various tree nodes
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.chem.serialization import deserialize_action, serialize_action
from aizynthfinder.search.andor_trees import TreeNodeMixin
from aizynthfinder.search.retrostar.cost import MoleculeCost

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.chem.serialization import (
        MoleculeDeserializer,
        MoleculeSerializer,
    )
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import List, Optional, Sequence, Set, StrDict


class MoleculeNode(TreeNodeMixin):
    """
    An OR node representing a molecule

    :ivar cost: the cost of synthesizing the molecule
    :ivar expandable: if True, this node is part of the frontier
    :ivar mol: the molecule represented by the node
    :ivar in_stock: if True the molecule is in stock and hence should not be expanded
    :ivar parent: the parent of the node
    :ivar solved: if True the molecule is in stock or at least one child node is solved
    :ivar value: the current rn(m|T)

    :param mol: the molecule to be represented by the node
    :param config: the configuration of the search
    :param parent: the parent of the node, optional
    """

    def __init__(
        self,
        mol: TreeMolecule,
        config: Configuration,
        molecule_cost: MoleculeCost,
        parent: Optional[ReactionNode] = None,
    ) -> None:
        self.mol = mol
        self._config = config
        self.molecule_cost = molecule_cost
        self.cost = self.molecule_cost(mol)
        self.value = self.cost
        self.in_stock = mol in config.stock
        self.parent = parent

        self._children: List[ReactionNode] = []
        self.solved = self.in_stock
        # Makes it unexpandable if we have reached maximum depth
        self.expandable = self.mol.transform < self._config.search.max_transforms

        if self.in_stock:
            self.expandable = False
            self.value = 0

    @classmethod
    def create_root(
        cls, smiles: str, config: Configuration, molecule_cost: MoleculeCost
    ) -> "MoleculeNode":
        """
        Create a root node for a tree using a SMILES.

        :param smiles: the SMILES representation of the root state
        :param config: settings of the tree search algorithm
        :return: the created node
        """
        mol = TreeMolecule(parent=None, transform=0, smiles=smiles)
        return MoleculeNode(mol=mol, config=config, molecule_cost=molecule_cost)

    @classmethod
    def from_dict(
        cls,
        dict_: StrDict,
        config: Configuration,
        molecules: MoleculeDeserializer,
        molecule_cost: MoleculeCost,
        parent: Optional[ReactionNode] = None,
    ) -> "MoleculeNode":
        """
        Create a new node from a dictionary, i.e. deserialization

        :param dict_: the serialized node
        :param config: settings of the tree search algorithm
        :param molecules: the deserialized molecules
        :param parent: the parent node
        :return: a deserialized node
        """
        mol = molecules.get_tree_molecules([dict_["mol"]])[0]
        node = MoleculeNode(mol, config, molecule_cost, parent)
        node.molecule_cost = molecule_cost
        for attr in ["cost", "expandable", "value"]:
            setattr(node, attr, dict_[attr])
        node.children = [
            ReactionNode.from_dict(
                child, config, molecules, node.molecule_cost, parent=node
            )
            for child in dict_["children"]
        ]
        return node

    @property  # type: ignore
    def children(self) -> List[ReactionNode]:  # type: ignore
        """Gives the reaction children nodes"""
        return self._children

    @children.setter
    def children(self, value: List[ReactionNode]) -> None:
        self._children = value

    @property
    def target_value(self) -> float:
        """
        The V_t(m|T) value,
        the current cost of the tree containing this node

        :return: the value
        """
        if self.parent:
            return self.parent.target_value
        return self.value

    @property
    def prop(self) -> StrDict:
        return {"solved": self.solved, "mol": self.mol}

    def add_stub(self, cost: float, reaction: RetroReaction) -> Sequence[MoleculeNode]:
        """
        Add a stub / sub-tree to this node

        :param cost: the cost of the reaction
        :param reaction: the reaction creating the stub
        :return: list of all newly added molecular nodes
        """
        reactants = reaction.reactants[reaction.index]
        if not reactants:
            return []

        ancestors = self.ancestors()
        for mol in reactants:
            if mol in ancestors:
                return []

        rxn_node = ReactionNode.create_stub(
            cost=cost,
            reaction=reaction,
            parent=self,
            config=self._config,
        )
        self._children.append(rxn_node)

        return rxn_node.children

    def ancestors(self) -> Set[TreeMolecule]:
        """
        Return the ancestors of this node

        :return: the ancestors
        :rtype: set
        """
        if not self.parent:
            return {self.mol}

        ancestors = self.parent.parent.ancestors()
        ancestors.add(self.mol)
        return ancestors

    def close(self) -> float:
        """
        Updates the values of this node after expanding it.

        :return: the delta V value
        :rtype: float
        """
        self.solved = any(child.solved for child in self.children)
        if self.children:
            new_value = np.min([child.value for child in self.children])
        else:
            new_value = np.inf

        v_delta = new_value - self.value
        self.value = new_value

        self.expandable = False
        return v_delta

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :return: the serialized node
        """
        dict_ = {attr: getattr(self, attr) for attr in ["cost", "expandable", "value"]}
        dict_["mol"] = molecule_store[self.mol]
        dict_["children"] = [child.serialize(molecule_store) for child in self.children]
        return dict_

    def update(self, solved: bool) -> None:
        """
        Update the node as part of the update algorithm,
        calling the `update()` method of its parent if available.

        :param solved: if the child node was solved
        """
        new_value = np.min([child.value for child in self.children])
        new_solv = self.solved or solved
        updated = (self.value != new_value) or (self.solved != new_solv)

        v_delta = new_value - self.value
        self.value = new_value
        self.solved = new_solv

        if updated and self.parent:
            self.parent.update(v_delta, from_mol=self.mol)


class ReactionNode(TreeNodeMixin):
    """
    An AND node representing a reaction

    :ivar cost: the cost of the reaction
    :ivar parent: the parent of the node
    :ivar reaction: the reaction represented by the node
    :ivar solved: if True all children nodes are solved
    :ivar target_value: the V(m|T) for the children, the current cost
    :ivar value: the current rn(r|T)

    :param cost: the cost of the reaction
    :param reaction: the reaction to be represented by the node
    :param parent: the parent of the node
    """

    def __init__(
        self, cost: float, reaction: RetroReaction, parent: MoleculeNode
    ) -> None:
        self.parent = parent
        self.cost = cost
        self.reaction = reaction

        self._children: List[MoleculeNode] = []
        self.solved = False
        # rn(R|T)
        self.value = self.cost
        # V(R|T) = V(m|T) for m in children
        self.target_value = self.parent.target_value - self.parent.value + self.value

    @classmethod
    def create_stub(
        cls,
        cost: float,
        reaction: RetroReaction,
        parent: MoleculeNode,
        config: Configuration,
    ) -> ReactionNode:
        """
        Create a ReactionNode and creates all the MoleculeNode objects
        that are the children of the node.

        :param cost: the cost of the reaction
        :param reaction: the reaction to be represented by the node
        :param parent: the parent of the node
        :param config: the configuration of the search tree
        """
        node = cls(cost, reaction, parent)
        reactants = reaction.reactants[reaction.index]
        node.children = [
            MoleculeNode(
                mol=mol, config=config, molecule_cost=parent.molecule_cost, parent=node
            )
            for mol in reactants
        ]
        node.solved = all(child.solved for child in node.children)
        # rn(R|T)
        node.value = node.cost + sum(child.value for child in node.children)
        # V(R|T) = V(m|T) for m in children
        node.target_value = node.parent.target_value - node.parent.value + node.value
        return node

    @classmethod
    def from_dict(
        cls,
        dict_: StrDict,
        config: Configuration,
        molecules: MoleculeDeserializer,
        molecule_cost: MoleculeCost,
        parent: MoleculeNode,
    ) -> ReactionNode:
        """
        Create a new node from a dictionary, i.e. deserialization

        :param dict_: the serialized node
        :param config: the configuration of the tree search
        :param molecules: the deserialized molecules
        :param parent: the parent node
        :return: a deserialized node
        """
        reaction = deserialize_action(dict_["reaction"], molecules)
        node = cls(0, reaction, parent)
        for attr in ["cost", "value", "target_value"]:
            setattr(node, attr, dict_[attr])
        node.children = [
            MoleculeNode.from_dict(child, config, molecules, molecule_cost, parent=node)
            for child in dict_["children"]
        ]
        node.solved = all(child.solved for child in node.children)
        return node

    @property  # type: ignore
    def children(self) -> List[MoleculeNode]:  # type: ignore
        """Gives the molecule children nodes"""
        return self._children

    @children.setter
    def children(self, value: List[MoleculeNode]) -> None:
        self._children = value

    @property
    def prop(self) -> StrDict:
        return {"solved": self.solved, "reaction": self.reaction}

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :return: the serialized node
        """
        dict_ = {
            attr: getattr(self, attr) for attr in ["cost", "value", "target_value"]
        }
        dict_["reaction"] = serialize_action(self.reaction, molecule_store)
        dict_["children"] = [child.serialize(molecule_store) for child in self.children]
        return dict_

    def update(self, value: float, from_mol: Optional[TreeMolecule] = None) -> None:
        """
        Update the node as part of the update algorithm,
        calling the `update()` method of its parent

        :param value: the delta V value
        :param from_mol: the molecule being expanded, used for excluding propagation
        """
        self.value += value
        self.target_value += value
        self.solved = all(node.solved for node in self.children)

        if value != 0:
            self._propagate(value, exclude=from_mol)

        self.parent.update(self.solved)

    def _propagate(self, value: float, exclude: Optional[TreeMolecule] = None) -> None:
        if not exclude:
            self.target_value += value

        for child in self.children:
            if exclude is None or child.mol is not exclude:
                for grandchild in child.children:
                    grandchild._propagate(value)  # pylint: disable=protected-access
