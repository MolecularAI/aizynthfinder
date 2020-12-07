"""
Helper routines and class for the `aizynthfinder.analysis` module.
To avoid clutter in that module, larger utility algorithms are placed herein.
"""
import abc
from collections import defaultdict

import networkx as nx

from aizynthfinder.chem import (
    Molecule,
    UniqueMolecule,
    FixedRetroReaction,
    hash_reactions,
)
from aizynthfinder.utils.image import make_visjs_page


class _ReactionTreeLoader(abc.ABC):
    """ Base class for classes that creates a reaction tree object
    """

    def __init__(self, *args, **kwargs):
        # To avoid circular imports
        from aizynthfinder.analysis import ReactionTree  # noqa

        self.tree = ReactionTree()
        self._load(*args, **kwargs)

        self.tree.is_solved = all(
            self.tree.in_stock(node) for node in self.tree.leafs()
        )
        _RepeatingPatternIdentifier.find(self.tree)

    def _add_node(self, node, depth=0, transform=0, in_stock=False, hide=False):
        attributes = {
            "hide": hide,
            "depth": depth,
        }
        if isinstance(node, Molecule):
            attributes.update({"transform": transform, "in_stock": in_stock})
        self.tree.graph.add_node(node, **attributes)
        if not self.tree.root:
            self.tree.root = node

    @abc.abstractmethod
    def _load(self, *args, **kwargs):
        pass


class ReactionTreeFromDict(_ReactionTreeLoader):
    """
    Creates a reaction tree object from a dictionary

    :param tree_dict: the dictionary representation
    :type tree_dict: dict
    """

    def _load(self, tree_dict):
        self._parse_tree_dict(tree_dict)

    def _parse_tree_dict(self, tree_dict, ncalls=0):
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
        reaction_node.reactants = reactant_nodes

        return product_node


class ReactionTreeFromMcts(_ReactionTreeLoader):
    """
    Creates a reaction tree object from MCTS nodes and reaction objects

    :param actions: the reactions forming the route
    :type actions: list of RetroReaction
    :param nodes: the MCTS nodes forming the route
    :type nodes: list of Node
    """

    def _load(self, actions, nodes):
        self._unique_mols = {}

        root_mol = nodes[0].state.mols[0]
        self._unique_mols[id(root_mol)] = root_mol.make_unique()
        self._add_node(
            self._unique_mols[id(root_mol)], in_stock=nodes[0].state.is_solved,
        )

        for child, action in zip(nodes[1:], actions):
            self._add_bipartite(child, action)

    def _add_bipartite(self, child, action):

        reaction_obj = FixedRetroReaction(
            self._unique_mol(action.mol), smiles=action.smiles, metadata=action.metadata
        )
        self._add_node(reaction_obj, depth=action.mol.transform + 1)
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
        reaction_obj.reactants = reactant_nodes

    def _unique_mol(self, molecule):
        id_ = id(molecule)
        if id_ not in self._unique_mols:
            self._unique_mols[id_] = molecule.make_unique()
        return self._unique_mols[id_]


class _RepeatingPatternIdentifier:
    """
    Encapsulation of algorithm to identify repeating patterns of reactions and mark them as hidden.

    A unit of the repetition is the hash of two consecutive reactions,
    where the first unit should be the first two reactions of the route.

    This is for hiding repeating patterns of e.g. protection followed by deprotection,
    which is a common behaviour for the tree search when it fails to solve a route.
    """

    @staticmethod
    def find(reaction_tree):
        """
        Find the repeating patterns and mark the nodes

        :param reaction_tree: the reaction tree to process
        :type reaction_tree: ReactionTree
        """
        for node in reaction_tree.reactions():
            # We are only interesting of starting at the very first reaction
            if any(reaction_tree.graph[mol] for mol in node.reactants[0]):
                continue
            actions = _RepeatingPatternIdentifier._list_reactions(reaction_tree, node)
            if len(actions) < 5:
                continue

            hashes = [
                hash_reactions([rxn1, rxn2], sort=False)
                for rxn1, rxn2 in zip(actions[:-1:2], actions[1::2])
            ]
            for idx, (hash1, hash2) in enumerate(zip(hashes[:-1], hashes[1:])):
                if hash1 == hash2:
                    _RepeatingPatternIdentifier._hide_reaction(
                        reaction_tree, actions[idx * 2]
                    )
                    _RepeatingPatternIdentifier._hide_reaction(
                        reaction_tree, actions[idx * 2 + 1]
                    )
                    reaction_tree.has_repeating_patterns = True
                # The else-clause prevents removing repeating patterns in the middle of a route
                else:
                    break

    @staticmethod
    def _hide_reaction(reaction_tree, reaction_node):
        reaction_tree.graph.nodes[reaction_node]["hide"] = True
        for reactants in reaction_node.reactants[0]:
            reaction_tree.graph.nodes[reactants]["hide"] = True

    @staticmethod
    def _list_reactions(reaction_tree, reaction_node):
        """ List all reaction nodes from the given one to the last
        """
        reactions = [reaction_node]
        curr_rxn = reaction_node
        product = reaction_node.mol
        while product is not reaction_tree.root:
            curr_rxn = next(reaction_tree.graph.predecessors(product))
            product = curr_rxn.mol
            reactions.append(curr_rxn)
        return reactions


class CombinedReactionTrees:
    """
    Encapsulation of an algorithm that combines several reaction trees into a
    larger bipartite graph with all reactions and molecules.

    The reactions at a specific level of the reaction trees are grouped based
    on the reaction smiles.

    :params reactions_trees: the list of reaction trees to combine
    :type reaction_trees: list of ReactionTree
    """

    def __init__(self, reaction_trees):
        self.graph = nx.DiGraph()
        first_rt = reaction_trees[0]
        # This is to avoid circular imports
        self._reaction_tree_class = first_rt.__class__
        self.root = first_rt.root
        self.graph.add_node(self.root, in_stock=first_rt.in_stock(self.root))
        rt_node_spec = [(rt.root, rt.graph) for rt in reaction_trees]
        self._add_reaction_trees_to_node(self.root, rt_node_spec)

    def to_dict(self):
        """
        Returns the graph as a dictionary in a pre-defined format.

        :return: the combined reaction trees
        :rtype: dict
        """
        rt = self._reaction_tree_class()
        rt.root = self.root
        rt.graph = self.graph
        return rt.to_dict()

    def to_visjs_page(
        self, filename, in_stock_colors={True: "green", False: "orange"},
    ):
        """
        Create a visualization of the combined reaction tree using the vis.js network library.

        The HTML page and all the images will be put into a tar-ball.

        :param filename: the name of the tarball
        :type filename: str
        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :type in_stock_colors: dict, optional
        """
        molecules = [node for node in self.graph if isinstance(node, Molecule)]
        reactions = [node for node in self.graph if not isinstance(node, Molecule)]
        frame_colors = [
            in_stock_colors[self.graph.nodes[node].get("in_stock", False)]
            for node in molecules
        ]
        make_visjs_page(filename, molecules, reactions, self.graph.edges, frame_colors)

    def _add_reaction_trees_to_node(self, base_node, rt_node_spec):

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
    def _find_other_children(child, group):
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
