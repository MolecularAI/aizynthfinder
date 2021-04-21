"""
Module containing helper classes to compute the distance between to reaction trees using the APTED method
Since APTED is based on ordered trees and the reaction trees are unordered, plenty of
heuristics are implemented to deal with this.
"""
from __future__ import annotations
import random
import itertools
import math
from enum import Enum
from operator import itemgetter
from typing import TYPE_CHECKING

import numpy as np
from apted import Config as BaseAptedConfig
from apted import APTED as Apted
from scipy.spatial.distance import jaccard as jaccard_dist

from aizynthfinder.chem import Molecule
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        List,
        Union,
        Iterable,
        Tuple,
        Dict,
        StrDict,
    )
    from aizynthfinder.analysis import ReactionTree
    from aizynthfinder.chem import Reaction

    _FloatIterator = Iterable[float]


class TreeContent(str, Enum):
    """Possibilities for distance calculations on reaction trees"""

    MOLECULES = "molecules"
    REACTIONS = "reactions"
    BOTH = "both"


class AptedConfig(BaseAptedConfig):
    """
    This is a helper class for the tree edit distance
    calculation. It defines how the substitution
    cost is calculated and how to obtain children nodes.

    :param randomize: if True, the children will be shuffled
    :param sort_children: if True, the children will be sorted
    """

    def __init__(self, randomize: bool = False, sort_children: bool = False) -> None:
        super().__init__()
        self._randomize = randomize
        self._sort_children = sort_children

    def rename(self, node1: StrDict, node2: StrDict) -> float:
        if node1["type"] != node2["type"]:
            return 1

        fp1 = node1["fingerprint"]
        fp2 = node2["fingerprint"]
        return jaccard_dist(fp1, fp2)

    def children(self, node: StrDict) -> List[StrDict]:
        if self._sort_children:
            return sorted(node["children"], key=itemgetter("sort_key"))
        if not self._randomize:
            return node["children"]
        children = list(node["children"])
        random.shuffle(children)
        return children


class ReactionTreeWrapper:
    """
    Wrapper or a reaction tree that can calculate distances between
    trees.

    :param reaction_tree: the reaction tree to wrap
    :param content: the content of the route to consider in the distance calculation
    :param exhaustive_limit: if the number of possible ordered trees are below this limit create them all
    """

    _index_permutations = {
        n: list(itertools.permutations(range(n), n)) for n in range(1, 8)
    }

    def __init__(
        self,
        reaction_tree: ReactionTree,
        content: Union[str, TreeContent] = TreeContent.MOLECULES,
        exhaustive_limit: int = 20,
    ) -> None:
        self._logger = logger()
        # Will convert string input automatically
        self._content = TreeContent(content)
        self._graph = reaction_tree.graph
        self._root = self._make_root(reaction_tree)

        self._trees = []
        self._tree_count, self._node_index_list = self._inspect_tree()
        self._enumeration = self._tree_count <= exhaustive_limit

        if not self._root:
            return

        if self._enumeration:
            self._create_all_trees()
        else:
            self._trees.append(self._create_tree_recursively(self._root))

    @property
    def info(self) -> StrDict:
        """Return a dictionary with internal information about the wrapper"""
        return {
            "content": self._content,
            "tree count": self._tree_count,
            "enumeration": self._enumeration,
            "root": self._root,
        }

    @property
    def first_tree(self) -> StrDict:
        """Return the first created ordered tree"""
        return self._trees[0]

    @property
    def trees(self) -> List[StrDict]:
        """Return a list of all created ordered trees"""
        return self._trees

    def distance_iter(
        self, other: "ReactionTreeWrapper", exhaustive_limit: int = 20
    ) -> _FloatIterator:
        """
        Iterate over all distances computed between this and another tree

        There are three possible enumeration of distances possible dependent
        on the number of possible ordered trees for the two routes that are compared

        * If the product of the number of possible ordered trees for both routes are
          below `exhaustive_limit` compute the distance between all pair of trees
        * If both self and other has been fully enumerated (i.e. all ordered trees has been created)
          compute the distances between all trees of the route with the most ordered trees and
          the first tree of the other route
        * Compute `exhaustive_limit` number of distances by shuffling the child order for
          each of the routes.

        The rules are applied top-to-bottom.

        :param other: another tree to calculate distance to
        :param exhaustive_limit: used to determine what type of enumeration to do
        :yield: the next computed distance between self and other
        """
        if self._tree_count * other.info["tree count"] < exhaustive_limit:
            yield from self._distance_iter_exhaustive(other)
        elif self._enumeration or other.info["enumeration"]:
            yield from self._distance_iter_semi_exhaustive(other)
        else:
            yield from self._distance_iter_random(other, exhaustive_limit)

    def distance_to(
        self, other: "ReactionTreeWrapper", exhaustive_limit: int = 20
    ) -> float:
        """
        Calculate the minimum distance from this route to another route

        Enumerate the distances using `distance_iter`.

        :param other: another tree to calculate distance to
        :param exhaustive_limit: used to determine what type of enumeration to do
        :return: the minimum distance
        """
        min_dist = 1e6
        min_iter = -1
        for iteration, distance in enumerate(
            self.distance_iter(other, exhaustive_limit)
        ):
            if distance < min_dist:
                min_iter = iteration
                min_dist = distance
        self._logger.debug(f"Found minimum after {min_iter} iterations")
        return min_dist

    def distance_to_with_sorting(self, other: "ReactionTreeWrapper") -> float:
        """
        Compute the distance to another tree, by simpling sorting the children
        of both trees. This is not guaranteed to return the minimum distance.

        :param other: another tree to calculate distance to
        :return: the distance
        """
        config = AptedConfig(sort_children=True)
        return Apted(self.first_tree, other.first_tree, config).compute_edit_distance()

    def _compute_fingerprint(self, node: Union[Molecule, Reaction]) -> np.ndarray:
        if isinstance(node, Molecule):
            return node.fingerprint(radius=2).astype(int)

        # Difference fingerprint for reactions
        product = next(self._graph.predecessors(node))
        fp = product.fingerprint(radius=2).copy()
        for reactant in self._graph.successors(node):
            fp -= reactant.fingerprint(radius=2)
        return fp.astype(int)

    def _create_all_trees(self) -> None:
        if not self._root:
            return

        self._trees = []
        # Iterate over all possible combinations of child order
        for order_list in itertools.product(*self._node_index_list):
            order_dict = {idict["node"]: idict["child_order"] for idict in order_list}
            self._trees.append(self._create_tree_recursively(self._root, order_dict))

    def _create_tree_recursively(
        self,
        node: Union[Molecule, Reaction],
        order_dict: Dict[Union[Molecule, Reaction], List[int]] = None,
    ) -> StrDict:
        fp = self._compute_fingerprint(node)
        dict_tree: StrDict = {
            "type": node.__class__.__name__,
            "smiles": node.smiles,
            "fingerprint": fp,
            "sort_key": "".join(f"{digit}" for digit in fp),
            "children": [],
        }
        for child in self._iter_children(node, order_dict):
            child_tree = self._create_tree_recursively(child, order_dict)
            dict_tree["children"].append(child_tree)
        return dict_tree

    def _distance_iter_exhaustive(self, other: "ReactionTreeWrapper") -> _FloatIterator:
        self._logger.debug(
            f"APTED: Exhaustive search. {len(self.trees)} {len(other.trees)}"
        )
        config = AptedConfig(randomize=False)
        for tree1, tree2 in itertools.product(self.trees, other.trees):
            yield Apted(tree1, tree2, config).compute_edit_distance()

    def _distance_iter_random(
        self, other: "ReactionTreeWrapper", ntimes: int
    ) -> _FloatIterator:
        self._logger.debug(
            f"APTED: Heuristic search. {len(self.trees)} {len(other.trees)}"
        )
        config = AptedConfig(randomize=False)
        yield Apted(self.first_tree, other.first_tree, config).compute_edit_distance()

        config = AptedConfig(randomize=True)
        for _ in range(ntimes):
            yield Apted(
                self.first_tree, other.first_tree, config
            ).compute_edit_distance()

    def _distance_iter_semi_exhaustive(
        self, other: "ReactionTreeWrapper"
    ) -> _FloatIterator:
        self._logger.debug(
            f"APTED: Semi-exhaustive search. {len(self.trees)} {len(other.trees)}"
        )
        if len(self.trees) < len(other.trees):
            first_wrapper = self
            second_wrapper = other
        else:
            first_wrapper = other
            second_wrapper = self

        config = AptedConfig(randomize=False)
        for tree1 in first_wrapper.trees:
            yield Apted(
                tree1, second_wrapper.first_tree, config
            ).compute_edit_distance()

    def _inspect_tree(self) -> Tuple[int, List[List[StrDict]]]:
        """
        Find the number of children for each node in the tree, which
        will be used to compute the number of possible combinations of child orders

        Also accumulate the possible child orders for the nodes.
        """
        permutations = []
        node_index_list = []
        for node in self._graph.nodes:
            # fmt: off
            if (
                isinstance(node, Molecule) and self._content is TreeContent.REACTIONS
            ) or (
                not isinstance(node, Molecule) and self._content is TreeContent.MOLECULES
            ):
                continue
            # fmt: on

            nchildren = len(list(self._iter_children(node)))
            permutations.append(math.factorial(nchildren))

            if nchildren > 0:
                node_index_list.append(
                    [
                        {"node": node, "child_order": idx}
                        for idx in self._index_permutations[nchildren]
                    ]
                )
        if not permutations:
            return 0, []
        return int(np.prod(permutations)), node_index_list

    def _iter_children(
        self,
        node: Union[Molecule, Reaction],
        order_dict: Dict[Union[Molecule, Reaction], List[int]] = None,
    ) -> Iterable[Union[Molecule, Reaction]]:
        def _generator(node_, lookup_node):
            if order_dict is None:
                for child in self._graph.successors(node_):
                    yield child
            else:
                children = list(self._graph.successors(node_))
                if children:
                    for child_idx in order_dict.get(lookup_node, []):
                        yield children[child_idx]

        if self._content is TreeContent.BOTH:
            yield from _generator(node, node)
        else:
            for succ in self._graph.successors(node):
                yield from _generator(succ, node)

    def _make_root(
        self, reaction_tree: ReactionTree
    ) -> Union[Molecule, Reaction, None]:
        if self._content is TreeContent.REACTIONS:
            try:
                return next(self._graph.successors(reaction_tree.root))
            except StopIteration:
                return None
        else:
            return reaction_tree.root
