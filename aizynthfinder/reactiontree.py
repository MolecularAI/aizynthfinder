""" Module containing the implementation of a reaction tree or route and factory classes to make such trees """

from __future__ import annotations

import abc
import hashlib
import json
import operator
from typing import TYPE_CHECKING

import networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_tree
from rxnutils.routes.comparison import simple_route_similarity
from rxnutils.routes.image import RouteImageFactory
from rxnutils.routes.readers import read_aizynthfinder_dict

from aizynthfinder.chem import (
    FixedRetroReaction,
    Molecule,
    UniqueMolecule,
    none_molecule,
)

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        FrameColors,
        Iterable,
        List,
        Optional,
        PilImage,
        StrDict,
        Union,
    )


class ReactionTree:
    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either FixedRetroReaction or UniqueMolecule objects.

    The reaction tree is initialized at instantiation and is not supposed to
    be updated.

    :ivar graph: the bipartite graph
    :ivar is_solved: if all of the leaf nodes are in stock
    :ivar root: the root of the tree
    :ivar created_at_iteration: iteration the reaction tree was created
    """

    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.root = none_molecule()
        self.is_solved: bool = False
        self.created_at_iteration: Optional[int] = None

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

    @property
    def metadata(self) -> StrDict:
        """Return a dicitionary with route metadata"""
        return {
            "created_at_iteration": self.created_at_iteration,
            "is_solved": self.is_solved,
        }

    def child_reactions(self, reaction: FixedRetroReaction) -> List[FixedRetroReaction]:
        """
        Return the child reactions of a reaction node in the tree.

        :param reaction: the query node (reaction)
        :return: child reactions
        """
        child_molecule_nodes = self.graph.successors(reaction)
        reaction_nodes = []
        for molecule_node in child_molecule_nodes:
            reaction_nodes.extend(list(self.graph.successors(molecule_node)))
        return reaction_nodes

    def depth(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> int:
        """
        Return the depth of a node in the route

        :param node: the query node
        :return: the depth
        """
        return self.graph.nodes[node].get("depth", -1)

    def distance_to(self, other: "ReactionTree") -> float:
        """
        Calculate the distance to another reaction tree

        This is a simple distance based on a route similarity

        :param other: the reaction tree to compare to
        :return: the distance between the routes
        """
        route1 = read_aizynthfinder_dict(self.to_dict())
        route2 = read_aizynthfinder_dict(other.to_dict())
        return 1.0 - float(simple_route_similarity([route1, route2])[0, 1])

    def hash_key(self) -> str:
        """
        Calculates a hash code for the tree using the sha224 hash function recursively

        :return: the hash key
        """
        return self._hash_func(self.root)

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

    def parent_molecule(self, mol: UniqueMolecule) -> UniqueMolecule:
        """Returns the parent molecule within the reaction tree.
        :param mol: the query node (molecule)
        :return: the parent molecule
        """
        if mol is self.root:
            raise ValueError("Root molecule does not have any parent node.")

        parent_reaction = list(self.graph.predecessors(mol))[0]
        parent_molecule = list(self.graph.predecessors(parent_reaction))[0]
        return parent_molecule

    def reactions(self) -> Iterable[FixedRetroReaction]:
        """
        Generates the reaction nodes of the reaction tree

        :yield: the next reaction in the tree
        """
        for node in self.graph:
            if not isinstance(node, Molecule):
                yield node

    def subtrees(self) -> Iterable[ReactionTree]:
        """
        Generates the subtrees of this reaction tree a
        subtree is a reaction treee starting at a molecule node that has children.

        :yield: the next subtree
        """

        def create_subtree(root_node):
            subtree = ReactionTree()
            subtree.root = root_node
            subtree.graph = dfs_tree(self.graph, root_node)
            for node in subtree.graph:
                prop = dict(self.graph.nodes[node])
                prop["depth"] -= self.graph.nodes[root_node].get("depth", 0)
                if "transform" in prop:
                    prop["transform"] -= self.graph.nodes[root_node].get("transform", 0)
                subtree.graph.nodes[node].update(prop)
            subtree.is_solved = all(subtree.in_stock(node) for node in subtree.leafs())
            return subtree

        for node in self.molecules():
            if node is not self.root and self.graph[node]:
                yield create_subtree(node)

    def to_dict(self, include_metadata=False) -> StrDict:
        """
        Returns the reaction tree as a dictionary in a pre-defined format.
        :param include_metadata: if True include metadata
        :return: the reaction tree
        """
        return self._build_dict(self.root, include_metadata=include_metadata)

    def to_image(
        self,
        in_stock_colors: Optional[FrameColors] = None,
        show_all: bool = True,
    ) -> PilImage:
        """
        Return a pictorial representation of the route

        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :param show_all: if True, also show nodes that are marked as hidden
        :return: the image of the route
        """
        factory = RouteImageFactory(
            self.to_dict(), in_stock_colors=in_stock_colors, show_all=show_all
        )
        return factory.image

    def to_json(self, include_metadata=False) -> str:
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        """
        return json.dumps(
            self.to_dict(include_metadata=include_metadata), sort_keys=False, indent=2
        )

    def _build_dict(
        self,
        node: Union[UniqueMolecule, FixedRetroReaction],
        dict_: Optional[StrDict] = None,
        include_metadata=False,
    ) -> StrDict:
        if dict_ is None:
            dict_ = {}

        if node is self.root and include_metadata:
            dict_["route_metadata"] = self.metadata

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

    def _hash_func(self, node: Union[FixedRetroReaction, UniqueMolecule]) -> str:
        if isinstance(node, UniqueMolecule):
            hash_ = hashlib.sha224(node.inchi_key.encode())
        else:
            hash_ = hashlib.sha224(node.hash_key().encode())
        child_hashes = sorted(
            self._hash_func(child) for child in self.graph.successors(node)
        )
        for child_hash in child_hashes:
            hash_.update(child_hash.encode())
        return hash_.hexdigest()


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
            metadata = dict(reaction.metadata)
            if ":" in reaction.mapped_reaction_smiles():
                metadata["mapped_reaction_smiles"] = reaction.mapped_reaction_smiles()
            self._unique_reactions[id_] = FixedRetroReaction(
                self._unique_mol(reaction.mol),
                smiles=reaction.smiles,
                metadata=metadata,
            )
        return self._unique_reactions[id_]


class ReactionTreeFromDict(ReactionTreeLoader):
    """Creates a reaction tree object from a dictionary"""

    def _load(self, tree_dict: StrDict) -> None:  # type: ignore
        if tree_dict.get("route_metadata"):
            self.tree.created_at_iteration = tree_dict["route_metadata"].get(
                "created_at_iteration"
            )
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
