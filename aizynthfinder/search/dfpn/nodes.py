""" Module containing a classes representation various tree nodes
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.search.andor_trees import TreeNodeMixin

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.search.dfpn import SearchTree
    from aizynthfinder.utils.type_utils import List, Optional, Sequence, Set, StrDict

BIG_INT = int(1e10)


class _SuperNode(TreeNodeMixin):
    def __init__(self) -> None:
        # pylint: disable=invalid-name
        self.pn = 1  # Proof-number
        self.dn = 1  # Disproof-number
        self.pn_threshold = BIG_INT
        self.dn_threshold = BIG_INT
        self._children: List["_SuperNode"] = []
        self.expandable = True

    @property  # type: ignore
    def children(self) -> List[ReactionNode]:  # type: ignore
        """Gives the reaction children nodes"""
        return self._children  # type: ignore

    @property
    def closed(self) -> bool:
        """Return if the node is proven or disproven"""
        return self.proven or self.disproven

    @property
    def proven(self) -> bool:
        """Return if the node is proven"""
        return self.pn == 0

    @property
    def disproven(self) -> bool:
        """Return if the node is disproven"""
        return self.dn == 0

    def explorable(self) -> bool:
        """Return if the node can be explored by the search algorithm"""
        return not (
            self.closed or self.pn > self.pn_threshold or self.dn > self.dn_threshold
        )

    def reset(self) -> None:
        """Reset the thresholds"""
        if self.closed or self.expandable:
            return
        for child in self._children:
            child.reset()
        self.update()
        self.pn_threshold = BIG_INT
        self.dn_threshold = BIG_INT

    def update(self) -> None:
        """Update the proof and disproof numbers"""
        raise NotImplementedError("Implement a child class")

    def _set_disproven(self) -> None:
        self.pn = BIG_INT
        self.dn = 0

    def _set_proven(self) -> None:
        self.pn = 0
        self.dn = BIG_INT


class MoleculeNode(_SuperNode):
    """
    An OR node representing a molecule

    :ivar expandable: if True, this node is part of the frontier
    :ivar mol: the molecule represented by the node
    :ivar in_stock: if True the molecule is in stock and hence should not be expanded
    :ivar parent: the parent of the node
    :ivar pn: the proof number
    :ivar dn: the disproof number
    :ivar pn_threshold: the threshold for proof number
    :ivar dn_threshold: the threshold for disproof number

    :param mol: the molecule to be represented by the node
    :param config: the configuration of the search
    :param parent: the parent of the node, optional
    """

    def __init__(
        self,
        mol: TreeMolecule,
        config: Configuration,
        owner: SearchTree,
        parent: Optional[ReactionNode] = None,
    ) -> None:
        super().__init__()

        self.mol = mol
        self._config = config
        self.in_stock = mol in config.stock
        self.parent = parent
        self._edge_costs: List[int] = []
        self.tree = owner

        # Makes it unexpandable if we have reached maximum depth
        self.expandable = self.mol.transform < self._config.search.max_transforms

        if self.in_stock:
            self.expandable = False
            self._set_proven()
        elif not self.expandable:
            self._set_disproven()

    @classmethod
    def create_root(
        cls, smiles: str, config: Configuration, owner: SearchTree
    ) -> "MoleculeNode":
        """
        Create a root node for a tree using a SMILES.

        :param smiles: the SMILES representation of the root state
        :param config: settings of the tree search algorithm
        :return: the created node
        """
        mol = TreeMolecule(parent=None, transform=0, smiles=smiles)
        return MoleculeNode(mol=mol, config=config, owner=owner)

    @property
    def prop(self) -> StrDict:
        return {"solved": self.proven, "mol": self.mol}

    def expand(self) -> None:
        """Expand the molecule by utilising an expansion policy"""
        self.expandable = False
        reactions, priors = self._config.expansion_policy([self.mol])
        self.tree.profiling["expansion_calls"] += 1

        if not reactions:
            self._set_disproven()
            return

        costs = -np.log(np.clip(priors, 1e-3, 1.0))
        reaction_costs = []
        reactions_to_expand = []
        for reaction, cost in zip(reactions, costs):
            try:
                _ = reaction.reactants
                self.tree.profiling["reactants_generations"] += 1
            except:  # pylint: disable=bare-except
                continue
            if not reaction.reactants:
                continue
            for idx, _ in enumerate(reaction.reactants):
                rxn_copy = reaction.copy(idx)
                reactions_to_expand.append(rxn_copy)
                reaction_costs.append(cost)

        for cost, rxn in zip(reaction_costs, reactions_to_expand):
            self._add_child(rxn, cost)

        if not self._children:
            self._set_disproven()

    def promising_child(self) -> Optional[ReactionNode]:
        """
        Find and return the most promising child for exploration
        Updates the thresholds on that child
        """
        min_indices = np.argsort(
            [
                edge_cost + child.pn if not child.closed else BIG_INT
                for edge_cost, child in zip(self._edge_costs, self._children)
            ]
        )
        best_child = self._children[min_indices[0]]
        if len(self._children) > 1 and not self._children[min_indices[1]].closed:
            s2_pn = self._children[min_indices[1]].pn
        else:
            s2_pn = BIG_INT

        best_child.pn_threshold = (
            min(self.pn_threshold, s2_pn + 2) - self._edge_costs[min_indices[0]]
        )
        best_child.dn_threshold = self.dn_threshold - self.dn + best_child.dn
        return best_child

    def update(self) -> None:
        """Update the proof and disproof numbers"""
        func = all if self.parent is None else any
        if func(child.proven for child in self._children):
            self._set_proven()
            return
        if all(child.disproven for child in self._children):
            self._set_disproven()
            return

        child_dns = [child.dn for child in self._children if not child.closed]
        if not child_dns:
            self._set_proven()
            return

        self.dn = sum(child_dns)
        if self.dn >= BIG_INT:
            self.pn = 0
        else:
            self.pn = min(
                edge_cost + child.pn
                for edge_cost, child in zip(self._edge_costs, self._children)
                if not child.closed
            )
        return

    def _add_child(self, reaction: RetroReaction, _: float) -> None:
        reactants = reaction.reactants[reaction.index]
        if not reactants:
            return

        ancestors = self._ancestors()
        for mol in reactants:
            if mol in ancestors:
                return

        rxn_node = ReactionNode(
            reaction=reaction, config=self._config, owner=self.tree, parent=self
        )
        self._children.append(rxn_node)
        self._edge_costs.append(1)

    def _ancestors(self) -> Set[TreeMolecule]:
        if not self.parent:
            return {self.mol}

        # pylint: disable=protected-access
        ancestors = self.parent.parent._ancestors()
        ancestors.add(self.mol)
        return ancestors


class ReactionNode(_SuperNode):
    """
    An AND node representing a reaction

    :ivar parent: the parent of the node
    :ivar reaction: the reaction represented by the node
    :ivar pn: the proof number
    :ivar dn: the disproof number
    :ivar pn_threshold: the threshold for proof number
    :ivar dn_threshold: the threshold for disproof number
    :ivar expandable: if the node is expandable

    :param reaction: the reaction to be represented by the node
    :param config: the configuration of the search
    :param parent: the parent of the node
    """

    def __init__(
        self,
        reaction: RetroReaction,
        config: Configuration,
        owner: SearchTree,
        parent: MoleculeNode,
    ) -> None:
        super().__init__()
        self._config = config
        self.parent = parent
        self.reaction = reaction
        self.tree = owner

    @property  # type: ignore
    def children(self) -> List[MoleculeNode]:  # type: ignore
        """Gives the molecule children nodes"""
        return self._children  # type: ignore

    @property
    def prop(self) -> StrDict:
        return {"solved": self.proven, "reaction": self.reaction}

    @property
    def proven(self) -> bool:
        """Return if the node is proven"""
        if self.expandable:
            return False
        if self.pn == 0:
            return True
        return all(child.proven for child in self._children)

    @property
    def disproven(self) -> bool:
        """Return if the node is disproven"""
        if self.expandable:
            return False
        if self.dn == 0:
            return True
        return any(child.disproven for child in self._children)

    def expand(self) -> None:
        """Expand the node by creating nodes for each reactant"""
        self.expandable = False
        reactants = self.reaction.reactants[self.reaction.index]
        self._children = [
            MoleculeNode(mol=mol, config=self._config, owner=self.tree, parent=self)
            for mol in reactants
        ]

    def promising_child(self) -> Optional[MoleculeNode]:
        """
        Find and return the most promising child for exploration
        Updates the thresholds on that child
        """
        min_indices = np.argsort(
            [child.dn if not child.closed else BIG_INT for child in self._children]
        )

        best_child = self._children[min_indices[0]]
        if len(self._children) > 1 and not self._children[min_indices[1]].closed:
            s2_dn = self._children[min_indices[1]].dn
        else:
            s2_dn = BIG_INT

        best_child.pn_threshold = self.pn_threshold - self.pn + best_child.pn
        best_child.dn_threshold = min(self.dn_threshold, s2_dn + 1)
        return best_child

    def update(self) -> None:
        """Update the proof and disproof numbers"""
        if all(child.proven for child in self._children):
            self._set_proven()
            return
        if any(child.disproven for child in self._children):
            self._set_disproven()
            return

        self.pn = sum(child.pn for child in self._children if not child.closed)
        self.dn = min(child.dn for child in self._children if not child.closed)
