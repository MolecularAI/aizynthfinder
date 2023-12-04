""" Module containing a class that holds the tree search
"""
from __future__ import annotations

import json
from typing import TYPE_CHECKING

from aizynthfinder.chem.serialization import MoleculeDeserializer, MoleculeSerializer
from aizynthfinder.search.andor_trees import AndOrSearchTreeBase, SplitAndOrTree
from aizynthfinder.search.breadth_first.nodes import MoleculeNode
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.reactiontree import ReactionTree
    from aizynthfinder.utils.type_utils import List, Optional, Sequence


class SearchTree(AndOrSearchTreeBase):
    """
    Encapsulation of the a breadth-first exhaustive search algorithm

    :ivar config: settings of the tree search algorithm
    :ivar root: the root node

    :param config: settings of the tree search algorithm
    :param root_smiles: the root will be set to a node representing this molecule, defaults to None
    """

    def __init__(
        self, config: Configuration, root_smiles: Optional[str] = None
    ) -> None:
        super().__init__(config, root_smiles)
        self._mol_nodes: List[MoleculeNode] = []
        self._added_mol_nodes: List[MoleculeNode] = []
        self._logger = logger()

        if root_smiles:
            self.root: Optional[MoleculeNode] = MoleculeNode.create_root(
                root_smiles, config
            )
            self._mol_nodes.append(self.root)
        else:
            self.root = None

        self._routes: List[ReactionTree] = []

        self.profiling = {
            "expansion_calls": 0,
            "reactants_generations": 0,
        }

    @classmethod
    def from_json(cls, filename: str, config: Configuration) -> SearchTree:
        """
        Create a new search tree by deserialization from a JSON file

        :param filename: the path to the JSON node
        :param config: the configuration of the search tree
        :return: a deserialized tree
        """

        def _find_mol_nodes(node):
            for child_ in node.children:
                tree._mol_nodes.append(child_)  # pylint: disable=protected-access
                for grandchild in child_.children:
                    _find_mol_nodes(grandchild)

        tree = cls(config)
        with open(filename, "r") as fileobj:
            dict_ = json.load(fileobj)
        mol_deser = MoleculeDeserializer(dict_["molecules"])
        tree.root = MoleculeNode.from_dict(dict_["tree"], config, mol_deser)
        tree._mol_nodes.append(tree.root)  # pylint: disable=protected-access
        for child in tree.root.children:
            _find_mol_nodes(child)
        return tree

    @property
    def mol_nodes(self) -> Sequence[MoleculeNode]:  # type: ignore
        """Return the molecule nodes of the tree"""
        return self._mol_nodes

    def one_iteration(self) -> bool:
        """
        Perform one iteration expansion.
        Expands all expandable molecule nodes in the tree, which should be
        on the same depth of the tree.

        :raises StopIteration: if the search should be pre-maturely terminated
        :return: if a solution was found
        :rtype: bool
        """
        if self.root is None:
            raise ValueError("Root is undefined. Cannot make an iteration")

        self._routes = []
        self._added_mol_nodes = []

        for next_node in self._mol_nodes:
            if next_node.expandable:
                self._expand(next_node)

        if not self._added_mol_nodes:
            self._logger.debug("No new nodes added in breadth-first iteration")
            raise StopIteration

        self._mol_nodes.extend(self._added_mol_nodes)
        solved = all(node.in_stock for node in self._mol_nodes if not node.children)
        return solved

    def routes(self) -> List[ReactionTree]:
        """
        Extracts and returns routes from the AND/OR tree

        :return: the routes
        """
        if self.root is None:
            return []
        if not self._routes:
            self._routes = SplitAndOrTree(self.root, self.config.stock).routes
        return self._routes

    def serialize(self, filename: str) -> None:
        """
        Seralize the search tree to a JSON file

        :param filename: the path to the JSON file
        :type filename: str
        """
        if self.root is None:
            raise ValueError("Cannot serialize tree as root is not defined")

        mol_ser = MoleculeSerializer()
        dict_ = {"tree": self.root.serialize(mol_ser), "molecules": mol_ser.store}
        with open(filename, "w") as fileobj:
            json.dump(dict_, fileobj, indent=2)

    def _expand(self, node: MoleculeNode) -> None:
        node.expandable = False
        reactions, _ = self.config.expansion_policy([node.mol])
        self.profiling["expansion_calls"] += 1

        if not reactions:
            return

        reactions_to_expand = []
        for reaction in reactions:
            try:
                self.profiling["reactants_generations"] += 1
                _ = reaction.reactants
            except:  # pylint: disable=bare-except
                continue
            if not reaction.reactants:
                continue
            for idx, _ in enumerate(reaction.reactants):
                rxn_copy = reaction.copy(idx)
                reactions_to_expand.append(rxn_copy)

        for rxn in reactions_to_expand:
            new_nodes = node.add_stub(rxn)
            self._added_mol_nodes.extend(new_nodes)
