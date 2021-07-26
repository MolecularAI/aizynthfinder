""" Module for base classes for AND/OR trees and some tree utility code """
from __future__ import annotations
import random
import abc
from typing import TYPE_CHECKING

import networkx as nx

from aizynthfinder.reactiontree import ReactionTree, ReactionTreeLoader
from aizynthfinder.chem import UniqueMolecule

if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.context.stock import Stock
    from aizynthfinder.chem import FixedRetroReaction
    from aizynthfinder.utils.type_utils import StrDict, Any, List, Union, Optional


class TreeNodeMixin:
    """ A mixin class for node in a tree """

    @property
    def prop(self) -> StrDict:
        """ Dictionary with publicly exposed properties """
        return {}

    @property
    def children(self) -> List["TreeNodeMixin"]:
        """ List of children nodes """
        return []


class AndOrSearchTreeBase(abc.ABC):
    """ A base class for a search tree based on an AND/OR structure """

    def __init__(self, config: Configuration, root_smiles: str = None) -> None:
        self.config = config
        self._root_smiles = root_smiles

    @property
    def mol_nodes(self) -> List[TreeNodeMixin]:
        """ Return the molecule nodes of the tree """
        return []

    @abc.abstractmethod
    def one_iteration(self) -> bool:  # pylint: disable=no-self-use
        """ Perform one iteration of the search """
        return False

    @abc.abstractmethod
    def routes(self) -> List[ReactionTree]:  # pylint: disable=no-self-use
        """ Return the routes of the tree """
        return []


class SplitAndOrTree:
    """
    Encapsulation of an algorithm to split an AND/OR tree into separate routes

    This is a modified algorithm of the one detailed in the CompRet paper:
    Shibukawa et al. (2020) J. Cheminf. 12, 52

    This implementation sets an upper-bound on the number of extracted routes
    to avoid combinatorial explosion.

    The routes are extracted on instantiation and the routes can be access from
    the `routes` attribute.

    :param root_node: the root of the AND/OR tree
    :param stock: the stock of the search
    :param max_routes: the maximum number of routes to extract
    """

    def __init__(
        self, root_node: TreeNodeMixin, stock: Stock, max_routes: int = 25000
    ) -> None:
        self._traces: List[_AndOrTrace] = []
        self._black_list: List[TreeNodeMixin] = []
        graph = _AndOrTrace(root_node)
        self._samples = {child: 0 for child in root_node.children}
        if root_node.children:
            self._sampling_cutoff = max_routes / len(root_node.children)
        else:
            self._sampling_cutoff = max_routes
        self._partition_search_tree(graph, root_node)
        self.routes = [
            ReactionTreeFromAndOrTrace(trace, stock).tree for trace in self._traces
        ]

    def _partition_search_tree(self, graph: _AndOrTrace, node: TreeNodeMixin) -> None:
        # fmt: off
        if self._sampling_cutoff and len(graph) > 1 and self._samples[graph.first_reaction] > self._sampling_cutoff:
            return
        # fmt: on

        if not node.children:
            self._traces.append(graph)

        children_to_search = [
            child for child in node.children if child not in self._black_list
        ]
        if not children_to_search:
            for child in node.children:
                self._black_list.remove(child)
            return

        graph_copy = graph.copy()

        child_node = self._select_child_node(children_to_search)
        graph.add_edge(node, child_node)
        for grandchild in child_node.children:
            graph.add_edge(child_node, grandchild)
        self._black_list.append(child_node)

        leaves = [node for node in graph if not graph[node] and node.children]
        if not leaves:
            self._traces.append(graph)
            self._samples[graph.first_reaction] += 1
        else:
            self._partition_search_tree(graph, leaves[0])

        self._partition_search_tree(graph_copy, node)

    def _select_child_node(self, children: List[TreeNodeMixin]) -> TreeNodeMixin:
        # This is what makes this algorithm different from what was
        # detailed in the CompRet paper
        if not self._sampling_cutoff:
            return children[0]

        solved_children = [
            child for child in children if child.prop.get("solved", False)
        ]
        if solved_children:
            return random.choice(solved_children)
        return random.choice(children)


class _AndOrTrace(nx.DiGraph):
    """Helper class for the SplitAndOrTree class."""

    def __init__(self, root: TreeNodeMixin = None) -> None:
        super().__init__()
        self.root = root
        self._first_reaction: Optional[TreeNodeMixin] = None
        self._reaction_tree: Optional[Any] = None
        if root:
            self.add_node(root)

    @property
    def first_reaction(self) -> TreeNodeMixin:
        """ Return the first reaction or raise an exception """
        assert self._first_reaction is not None
        return self._first_reaction

    def add_edge(
        self, u_of_edge: TreeNodeMixin, v_of_edge: TreeNodeMixin, **attr: Any
    ) -> None:
        if u_of_edge is self.root:
            if self._first_reaction is not None:
                raise ValueError(
                    "Re-defining trace. Trying to set first reaction twice."
                )
            self._first_reaction = v_of_edge
        super().add_edge(u_of_edge, v_of_edge, **attr)

    def copy(self, as_view: bool = False) -> "_AndOrTrace":
        other = super().copy(as_view)
        assert isinstance(other, _AndOrTrace)
        other.root = self.root
        other._first_reaction = self._first_reaction  # pylint: disable=protected-access
        return other


class ReactionTreeFromAndOrTrace(ReactionTreeLoader):
    """Creates a reaction tree object from an AND/OR Trace"""

    def _load(self, andor_trace: nx.DiGraph, stock: Stock) -> None:  # type: ignore
        """
        :param andor_trace: the trace graph
        :param stock: stock object
        """
        self._stock = stock
        self._trace_graph = andor_trace
        self._trace_root = self._find_root()

        self._add_node(
            self._unique_mol(self._trace_root.prop["mol"]),
            depth=0,
            transform=0,
            in_stock=self._trace_root.prop["mol"] in self._stock,
        )
        for node1, node2 in andor_trace.edges():
            rt_node1 = self._make_rt_node(node1)
            rt_node2 = self._make_rt_node(node2)
            self.tree.graph.add_edge(rt_node1, rt_node2)

    def _add_node_with_depth(
        self, node: Union[UniqueMolecule, FixedRetroReaction], base_node: TreeNodeMixin
    ) -> None:
        if node in self.tree.graph:
            return

        depth = nx.shortest_path_length(self._trace_graph, self._trace_root, base_node)
        if isinstance(node, UniqueMolecule):
            self._add_node(
                node, depth=depth, transform=depth // 2, in_stock=node in self._stock
            )
        else:
            self._add_node(node, depth=depth)

    def _find_root(self) -> TreeNodeMixin:
        for node, degree in self._trace_graph.in_degree():  # type: ignore
            if degree == 0:
                return node
        raise ValueError("Could not find root!")

    def _make_rt_node(
        self, node: TreeNodeMixin
    ) -> Union[UniqueMolecule, FixedRetroReaction]:
        if "mol" in node.prop:
            unique_obj = self._unique_mol(node.prop["mol"])
            self._add_node_with_depth(unique_obj, node)
            return unique_obj

        reaction_obj = self._unique_reaction(node.prop["reaction"])
        reaction_obj.reactants = (
            tuple(self._unique_mol(child.prop["mol"]) for child in node.children),
        )
        self._add_node_with_depth(reaction_obj, node)
        return reaction_obj
