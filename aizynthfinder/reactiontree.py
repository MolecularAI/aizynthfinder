""" Module containing the implementation of a reaction tree or route and factory classes to make such trees """
from __future__ import annotations
import operator
import json
import abc
from typing import TYPE_CHECKING

import networkx as nx
from route_distances.route_distances import route_distances_calculator

from aizynthfinder.chem import (
    Molecule,
    UniqueMolecule,
    FixedRetroReaction,
    none_molecule,
)
from aizynthfinder.utils.image import make_graphviz_image


if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        StrDict,
        PilImage,
        FrameColors,
        Union,
        Iterable,
        Any,
        Dict,
    )
    from aizynthfinder.chem import RetroReaction


class ReactionTree:
    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either FixedRetroReaction or UniqueMolecule objects.

    The reaction tree is initialized at instantiation and is not supposed to
    be updated.

    :ivar graph: the bipartite graph
    :ivar is_solved: if all of the leaf nodes are in stock
    :ivar root: the root of the tree
    """

    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.root = none_molecule()
        self.is_solved: bool = False

    @classmethod
    def from_dict(cls, tree_dict: StrDict) -> "ReactionTree":
        """
        Create a new ReactionTree by parsing a dictionary.

        This is supposed to be the opposite of ``to_dict``,
        but because that format loses information, the returned
        object is not a full copy as the stock will only contain
        the list of molecules marked as ``in_stock`` in the dictionary.

        The returned object should be sufficient to e.g. generate an image of the route.

        :param tree_dict: the dictionary representation
        :returns: the reaction tree
        """
        return ReactionTreeFromDict(tree_dict).tree

    def depth(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> int:
        """
        Return the depth of a node in the route

        :param node: the query node
        :return: the depth
        """
        return self.graph.nodes[node].get("depth", -1)

    def distance_to(self, other: "ReactionTree", content: str = "both") -> float:
        """
        Calculate the distance to another reaction tree

        This is a tree edit distance, with unit cost to
        insert and deleted nodes, and the Jaccard distance for substituting nodes

        :param other: the reaction tree to compare to
        :param content: determine what part of the tree to include in the calculation
        :return: the distance between the routes
        """
        calculator = route_distances_calculator("ted", content=content)
        distances = calculator([self.to_dict(), other.to_dict()])
        return distances[0, 1]

    def in_stock(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> bool:
        """
        Return if a node in the route is in stock

        Note that is a property set on creation and as such is not updated.

        :param node: the query node
        :return: if the molecule is in stock
        """
        return self.graph.nodes[node].get("in_stock", False)

    def is_branched(self) -> bool:
        """
        Returns if the route is branched

        i.e. checks if the maximum depth is not equal to the number
        of reactions.
        """
        nsteps = len(list(self.reactions()))
        max_depth = max(self.depth(leaf) for leaf in self.leafs())
        return nsteps != max_depth // 2

    def leafs(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecules nodes of the reaction tree that has no predecessors,
        i.e. molecules that has not been broken down

        :yield: the next leaf molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule) and not self.graph[node]:
                yield node

    def molecules(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecule nodes of the reaction tree

        :yield: the next molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule):
                yield node

    def reactions(self) -> Iterable[FixedRetroReaction]:
        """
        Generates the reaction nodes of the reaction tree

        :yield: the next reaction in the tree
        """
        for node in self.graph:
            if not isinstance(node, Molecule):
                yield node

    def to_dict(self) -> StrDict:
        """
        Returns the reaction tree as a dictionary in a pre-defined format.

        :return: the reaction tree
        """
        return self._build_dict(self.root)

    def to_image(
        self,
        in_stock_colors: FrameColors = None,
        show_all: bool = True,
    ) -> PilImage:
        """
        Return a pictorial representation of the route

        :raises ValueError: if image could not be produced
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :param show_all: if True, also show nodes that are marked as hidden
        :return: the image of the route
        """

        def show(node_):
            return not self.graph.nodes[node_].get("hide", False)

        in_stock_colors = in_stock_colors or {True: "green", False: "orange"}
        molecules = []
        frame_colors = []
        mols = sorted(
            self.molecules(),
            key=lambda mol: (
                self.depth(mol),
                mol.weight,
            ),
        )
        for node in mols:
            if not show_all and not show(node):
                continue
            molecules.append(node)
            frame_colors.append(in_stock_colors[self.in_stock(node)])

        reactions = [node for node in self.reactions() if show_all or show(node)]
        edges = [
            (node1, node2)
            for node1, node2 in self.graph.edges()
            if show_all or (show(node1) and show(node2))
        ]

        try:
            return make_graphviz_image(molecules, reactions, edges, frame_colors)
        except FileNotFoundError as err:
            raise ValueError(str(err))

    def to_json(self) -> str:
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        """
        return json.dumps(self.to_dict(), sort_keys=False, indent=2)

    def _build_dict(
        self, node: Union[UniqueMolecule, FixedRetroReaction], dict_: StrDict = None
    ) -> StrDict:
        if dict_ is None:
            dict_ = {}

        dict_["type"] = "mol" if isinstance(node, Molecule) else "reaction"
        dict_["hide"] = self.graph.nodes[node].get("hide", False)
        dict_["smiles"] = node.smiles
        if isinstance(node, UniqueMolecule):
            dict_["is_chemical"] = True
            dict_["in_stock"] = self.in_stock(node)
        elif isinstance(node, FixedRetroReaction):
            dict_["is_reaction"] = True
            dict_["metadata"] = dict(node.metadata)
        else:
            raise ValueError(
                f"This is an invalid reaction tree. Unknown node type {type(node)}"
            )

        dict_["children"] = []

        children = list(self.graph.successors(node))
        if isinstance(node, FixedRetroReaction):
            children.sort(key=operator.attrgetter("weight"))
        for child in children:
            child_dict = self._build_dict(child)
            dict_["children"].append(child_dict)

        if not dict_["children"]:
            del dict_["children"]
        return dict_


class ReactionTreeLoader(abc.ABC):
    """
    Base class for classes that creates a reaction tree object

    This class makes sure that node attributes are set after the
    graph is generated, and provides utility methods.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._unique_mols: Dict[int, UniqueMolecule] = {}
        self._unique_reactions: Dict[int, FixedRetroReaction] = {}
        self.tree = ReactionTree()
        self._load(*args, **kwargs)

        self.tree.is_solved = all(
            self.tree.in_stock(node) for node in self.tree.leafs()
        )

    def _add_node(
        self,
        node: Union[UniqueMolecule, FixedRetroReaction],
        depth: int = 0,
        transform: int = 0,
        in_stock: bool = False,
        hide: bool = False,
    ) -> None:
        attributes = {
            "hide": hide,
            "depth": depth,
        }
        if isinstance(node, Molecule):
            attributes.update({"transform": transform, "in_stock": in_stock})
            if not self.tree.root:
                self.tree.root = node
        self.tree.graph.add_node(node, **attributes)

    @abc.abstractmethod
    def _load(self, *args: Any, **kwargs: Any) -> None:
        pass

    def _unique_mol(self, molecule: Molecule) -> UniqueMolecule:
        id_ = id(molecule)
        if id_ not in self._unique_mols:
            self._unique_mols[id_] = molecule.make_unique()
        return self._unique_mols[id_]

    def _unique_reaction(self, reaction: RetroReaction) -> FixedRetroReaction:
        id_ = id(reaction)
        if id_ not in self._unique_reactions:
            self._unique_reactions[id_] = FixedRetroReaction(
                self._unique_mol(reaction.mol),
                smiles=reaction.smiles,
                metadata=reaction.metadata,
            )
        return self._unique_reactions[id_]


class ReactionTreeFromDict(ReactionTreeLoader):
    """Creates a reaction tree object from a dictionary"""

    def _load(self, tree_dict: StrDict) -> None:  # type: ignore
        self._parse_tree_dict(tree_dict)

    def _parse_tree_dict(self, tree_dict: StrDict, ncalls: int = 0) -> UniqueMolecule:
        product_node = UniqueMolecule(smiles=tree_dict["smiles"])
        self._add_node(
            product_node,
            depth=2 * ncalls,
            transform=ncalls,
            hide=tree_dict.get("hide", False),
            in_stock=tree_dict["in_stock"],
        )

        rxn_tree_dict = tree_dict.get("children", [])
        if not rxn_tree_dict:
            return product_node

        rxn_tree_dict = rxn_tree_dict[0]
        reaction_node = FixedRetroReaction(
            product_node,
            smiles=rxn_tree_dict["smiles"],
            metadata=rxn_tree_dict.get("metadata", {}),
        )
        self._add_node(
            reaction_node, depth=2 * ncalls + 1, hide=rxn_tree_dict.get("hide", False)
        )
        self.tree.graph.add_edge(product_node, reaction_node)

        reactant_nodes = []
        for reactant_tree in rxn_tree_dict.get("children", []):
            reactant_node = self._parse_tree_dict(reactant_tree, ncalls + 1)
            self.tree.graph.add_edge(reaction_node, reactant_node)
            reactant_nodes.append(reactant_node)
        reaction_node.reactants = (tuple(reactant_nodes),)

        return product_node


class ReactionTreeFromExpansion(ReactionTreeLoader):
    """
    Create a ReactionTree from a single reaction

    This is mainly intended as a convenience function for the expander interface
    """

    def _load(self, reaction: RetroReaction) -> None:  # type: ignore
        root = self._unique_mol(reaction.mol)
        self._add_node(root)

        rxn = self._unique_reaction(reaction)
        if hasattr(reaction, "smarts"):
            rxn.metadata["smarts"] = reaction.smarts  # type: ignore
        self._add_node(rxn)
        self.tree.graph.add_edge(root, rxn)

        reactant_nodes = []
        for reactant in reaction.reactants[0]:
            reactant_node = self._unique_mol(reactant)
            reactant_nodes.append(reactant_node)
            self._add_node(reactant_node)
            self.tree.graph.add_edge(rxn, reactant_node)
        rxn.reactants = (tuple(reactant_nodes),)
