""" Module containing a class that holds the tree search
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.andor_trees import AndOrSearchTreeBase, SplitAndOrTree
from aizynthfinder.search.dfpn.nodes import MoleculeNode, ReactionNode
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.search.andor_trees import TreeNodeMixin
    from aizynthfinder.utils.type_utils import List, Optional, Sequence, Union


class SearchTree(AndOrSearchTreeBase):
    """
    Encapsulation of the Depth-First Proof-Number (DFPN) search algorithm.

    This algorithm does not support:
        1. Filter policy
        2. Serialization and deserialization

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
        self._logger = logger()
        self._root_smiles = root_smiles
        if root_smiles:
            self.root: Optional[MoleculeNode] = MoleculeNode.create_root(
                root_smiles, config, self
            )
            self._mol_nodes.append(self.root)
        else:
            self.root = None

        self._routes: List[ReactionTree] = []
        self._frontier: Optional[Union[MoleculeNode, ReactionNode]] = None
        self._initiated = False

        self.profiling = {
            "expansion_calls": 0,
            "reactants_generations": 0,
        }

    @property
    def mol_nodes(self) -> Sequence[MoleculeNode]:  # type: ignore
        """Return the molecule nodes of the tree"""
        return self._mol_nodes

    def one_iteration(self) -> bool:
        """
        Perform one iteration of expansion.

        If possible expand the frontier node twice, i.e. expanding an OR
        node and then and AND node. If frontier not expandable step up in the
        tree and find a new frontier to expand.

        If a solution is found, mask that tree for exploration and start over.

        :raises StopIteration: if the search should be pre-maturely terminated
        :return: if a solution was found
        :rtype: bool
        """
        if not self._initiated:
            if self.root is None:
                raise ValueError("Root is undefined. Cannot make an iteration")

            self._routes = []
            self._frontier = self.root
        assert self.root is not None

        while True:
            # Expand frontier, should be OR node
            assert isinstance(self._frontier, MoleculeNode)
            expanded_or = self._search_step()
            expanded_and = False
            if self._frontier:
                # Expand frontier again, this time an AND node
                assert isinstance(self._frontier, ReactionNode)
                expanded_and = self._search_step()
            if (
                expanded_or
                or expanded_and
                or self._frontier is None
                or self._frontier is self.root
            ):
                break

        found_solution = any(child.proven for child in self.root.children)
        if self._frontier is self.root:
            self.root.reset()

        if self._frontier is None:
            raise StopIteration()

        return found_solution

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

    def _search_step(self) -> bool:
        assert self._frontier is not None
        expanded = False
        if self._frontier.expandable:
            self._frontier.expand()
            expanded = True
            if isinstance(self._frontier, ReactionNode):
                self._mol_nodes.extend(self._frontier.children)

        self._frontier.update()
        if not self._frontier.explorable():
            self._frontier = self._frontier.parent
            return False

        child = self._frontier.promising_child()
        if not child:
            self._frontier = self._frontier.parent
            return False

        self._frontier = child
        return expanded
