"""
Helper routines and class for the `aizynthfinder.analysis` package.
To avoid clutter in that package, larger utility algorithms are placed herein.
"""
from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from typing import TYPE_CHECKING

import networkx as nx

from aizynthfinder.chem import (
    Molecule,
    UniqueMolecule,
    FixedRetroReaction,
)
from aizynthfinder.utils.image import make_visjs_page
from aizynthfinder.reactiontree import ReactionTree

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Sequence,
        Tuple,
        StrDict,
        FrameColors,
    )


@dataclass
class RouteSelectionArguments:
    """
    Selection arguments for the tree analysis class

    If `return_all` is False, it will return at least `nmin` routes and if routes have the same
    score it will return them as well up to `nmax` routes.

    If `return_all` is True, it will return all solved routes if there is at least one is solved, otherwise
    the `nmin` and `nmax` will be used.
    """
    nmin: int = 5
    nmax: int = 25
    return_all: bool = False


class CombinedReactionTrees:
    """
    Encapsulation of an algorithm that combines several reaction trees into a
    larger bipartite graph with all reactions and molecules.

    The reactions at a specific level of the reaction trees are grouped based
    on the reaction smiles.

    :params reactions_trees: the list of reaction trees to combine
    """

    def __init__(self, reaction_trees: Sequence[ReactionTree]) -> None:
        self.graph = nx.DiGraph()
        first_rt = reaction_trees[0]
        self.root = first_rt.root

        self.graph.add_node(self.root, in_stock=first_rt.in_stock(self.root))
        rt_node_spec = [(rt.root, rt.graph) for rt in reaction_trees]
        self._add_reaction_trees_to_node(self.root, rt_node_spec)

    def to_dict(self) -> StrDict:
        """
        Returns the graph as a dictionary in a pre-defined format.

        :return: the combined reaction trees
        """
        rtree = ReactionTree()
        rtree.root = self.root
        rtree.graph = self.graph
        return rtree.to_dict()

    def to_visjs_page(
        self,
        filename: str,
        in_stock_colors: FrameColors = None,
    ) -> None:
        """
        Create a visualization of the combined reaction tree using the vis.js network library.

        The HTML page and all the images will be put into a tar-ball.

        :param filename: the name of the tarball
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        """
        in_stock_colors = in_stock_colors or {True: "green", False: "orange"}
        molecules = [node for node in self.graph if isinstance(node, Molecule)]
        reactions = [node for node in self.graph if not isinstance(node, Molecule)]
        frame_colors = [
            in_stock_colors[self.graph.nodes[node].get("in_stock", False)]
            for node in molecules
        ]
        make_visjs_page(filename, molecules, reactions, self.graph.edges, frame_colors)

    def _add_reaction_trees_to_node(
        self,
        base_node: UniqueMolecule,
        rt_node_spec: Sequence[Tuple[UniqueMolecule, nx.DiGraph]],
    ) -> None:

        reaction_groups = defaultdict(list)
        # Group the reactions from the nodes at this level based on the reaction smiles
        for node, graph in rt_node_spec:
            for reaction in graph[node]:
                reaction_groups[reaction.reaction_smiles()].append((graph, reaction))

        for group in reaction_groups.values():
            # Use the first RT in each group as the base
            first_graph, first_reaction = group[0]
            reaction_node = first_reaction.copy()
            self.graph.add_edge(base_node, reaction_node)

            for child in first_graph[first_reaction]:
                mol_node = child.make_unique()
                self.graph.add_node(
                    mol_node, in_stock=first_graph.nodes[child].get("in_stock", False)
                )
                self.graph.add_edge(reaction_node, mol_node)
                self._add_reaction_trees_to_node(
                    mol_node, self._find_other_children(child, group)
                )

    @staticmethod
    def _find_other_children(
        child: UniqueMolecule, group: Sequence[Tuple[nx.DiGraph, FixedRetroReaction]]
    ) -> Sequence[Tuple[UniqueMolecule, nx.DiGraph]]:
        children_spec = []
        for other_graph, other_reaction in group:
            found = False
            for other_child in other_graph[other_reaction]:
                if other_child.inchi_key == child.inchi_key:
                    children_spec.append((other_child, other_graph))
                    found = True
                    break
            if not found:
                raise ValueError("Could not find other child")
        return children_spec
