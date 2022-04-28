""" Module containing a classes representation various tree nodes
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.search.andor_trees import TreeNodeMixin
from aizynthfinder.chem.serialization import deserialize_action, serialize_action

if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem.serialization import (
        MoleculeDeserializer,
        MoleculeSerializer,
    )
    from aizynthfinder.utils.type_utils import (
        StrDict,
        Sequence,
        Set,
        List,
    )
    from aizynthfinder.chem import RetroReaction


class MoleculeNode(TreeNodeMixin):
    """
    An OR node representing a molecule

    :ivar expandable: if True, this node is part of the frontier
    :ivar mol: the molecule represented by the node
    :ivar in_stock: if True the molecule is in stock and hence should not be expanded
    :ivar parent: the parent of the node

    :param mol: the molecule to be represented by the node
    :param config: the configuration of the search
    :param parent: the parent of the node, optional
    """

    def __init__(
        self, mol: TreeMolecule, config: Configuration, parent: ReactionNode = None
    ) -> None:
        self.mol = mol
        self._config = config
        self.in_stock = mol in config.stock
        self.parent = parent

        self._children: List[ReactionNode] = []
        # Makes it unexpandable if we have reached maximum depth
        self.expandable = self.mol.transform <= self._config.max_transforms

        if self.in_stock:
            self.expandable = False

    @classmethod
    def create_root(cls, smiles: str, config: Configuration) -> "MoleculeNode":
        """
        Create a root node for a tree using a SMILES.

        :param smiles: the SMILES representation of the root state
        :param config: settings of the tree search algorithm
        :return: the created node
        """
        mol = TreeMolecule(parent=None, transform=0, smiles=smiles)
        return MoleculeNode(mol=mol, config=config)

    @classmethod
    def from_dict(
        cls,
        dict_: StrDict,
        config: Configuration,
        molecules: MoleculeDeserializer,
        parent: ReactionNode = None,
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
        node = MoleculeNode(mol, config, parent)
        node.expandable = dict_["expandable"]
        node.children = [
            ReactionNode.from_dict(child, config, molecules, parent=node)
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
    def prop(self) -> StrDict:
        return {"solved": self.in_stock, "mol": self.mol}

    def add_stub(self, reaction: RetroReaction) -> Sequence[MoleculeNode]:
        """
        Add a stub / sub-tree to this node

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
            reaction=reaction, parent=self, config=self._config
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

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :return: the serialized node
        """
        dict_: StrDict = {"expandable": self.expandable}
        dict_["mol"] = molecule_store[self.mol]
        dict_["children"] = [child.serialize(molecule_store) for child in self.children]
        return dict_


class ReactionNode(TreeNodeMixin):
    """
    An AND node representing a reaction

    :ivar parent: the parent of the node
    :ivar reaction: the reaction represented by the node

    :param cost: the cost of the reaction
    :param reaction: the reaction to be represented by the node
    :param parent: the parent of the node
    """

    def __init__(self, reaction: RetroReaction, parent: MoleculeNode) -> None:
        self.parent = parent
        self.reaction = reaction

        self._children: List[MoleculeNode] = []

    @classmethod
    def create_stub(
        cls,
        reaction: RetroReaction,
        parent: MoleculeNode,
        config: Configuration,
    ) -> ReactionNode:
        """
        Create a ReactionNode and creates all the MoleculeNode objects
        that are the children of the node.

        :param reaction: the reaction to be represented by the node
        :param parent: the parent of the node
        :param config: the configuration of the search tree
        """
        node = cls(reaction, parent)
        reactants = reaction.reactants[reaction.index]
        node.children = [
            MoleculeNode(mol=mol, config=config, parent=node) for mol in reactants
        ]
        return node

    @classmethod
    def from_dict(
        cls,
        dict_: StrDict,
        config: Configuration,
        molecules: MoleculeDeserializer,
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
        node = cls(reaction, parent)

        node.children = [
            MoleculeNode.from_dict(child, config, molecules, parent=node)
            for child in dict_["children"]
        ]
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
        return {"solved": False, "reaction": self.reaction}

    def serialize(self, molecule_store: MoleculeSerializer) -> StrDict:
        """
        Serialize the node object to a dictionary

        :param molecule_store: the serialized molecules
        :return: the serialized node
        """
        dict_ = {
            "reaction": serialize_action(self.reaction, molecule_store),
            "children": [child.serialize(molecule_store) for child in self.children],
        }
        return dict_
