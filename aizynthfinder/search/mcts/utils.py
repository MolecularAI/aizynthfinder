""" Module containing utility routines for MCTS. This is not part of public interface """
from __future__ import annotations
from typing import TYPE_CHECKING

from aizynthfinder.reactiontree import ReactionTreeLoader

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Tuple, List, Optional
    from aizynthfinder.search.mcts import MctsNode
    from aizynthfinder.chem import RetroReaction


class ReactionTreeFromSuperNode(ReactionTreeLoader):
    """
    Creates a reaction tree object from MCTS-like nodes and reaction objects
    """

    def _load(self, base_node: MctsNode) -> None:  # type: ignore
        actions, nodes = route_to_node(base_node)
        root_mol = nodes[0].state.mols[0]
        self._unique_mols[id(root_mol)] = root_mol.make_unique()
        self._add_node(
            self._unique_mols[id(root_mol)],
            in_stock=nodes[0].state.is_solved,
        )

        for child, action in zip(nodes[1:], actions):
            self._add_bipartite(child, action)

    def _add_bipartite(self, child: MctsNode, action: RetroReaction) -> None:

        reaction_obj = self._unique_reaction(action)
        self._add_node(reaction_obj, depth=2 * action.mol.transform + 1)
        self.tree.graph.add_edge(self._unique_mol(action.mol), reaction_obj)
        reactant_nodes = []
        for mol in child.state.mols:
            if mol.parent is action.mol:
                self._add_node(
                    self._unique_mol(mol),
                    depth=2 * mol.transform,
                    transform=mol.transform,
                    in_stock=mol in child.state.stock,
                )
                self.tree.graph.add_edge(reaction_obj, self._unique_mol(mol))
                reactant_nodes.append(self._unique_mol(mol))
        reaction_obj.reactants = (tuple(reactant_nodes),)


def route_to_node(
    from_node: MctsNode,
) -> Tuple[List[RetroReaction], List[MctsNode]]:
    """
    Return the route to a give node to the root.

    Will return both the actions taken to go between the nodes,
    and the nodes in the route themselves.

    :param from_node: the end of the route
    :return: the route
    """
    actions = []
    nodes = []
    current: Optional[MctsNode] = from_node

    while current is not None:
        parent = current.parent
        if parent is not None:
            action = parent[current]["action"]
            actions.append(action)
        nodes.append(current)
        current = parent
    return actions[::-1], nodes[::-1]
