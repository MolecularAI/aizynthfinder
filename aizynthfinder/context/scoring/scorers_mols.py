""" Module containing classes used to score the reaction routes.
These scores are only based on molecules in the route.
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from networkx import get_node_attributes
from rxnutils.chem.features.sc_score import SCScore

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.scoring.scorers_base import Scorer
from aizynthfinder.context.stock import StockException
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode

if TYPE_CHECKING:
    from aizynthfinder.chem import Molecule, UniqueMolecule
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.context.stock import Stock
    from aizynthfinder.utils.type_utils import (
        Iterable,
        List,
        Optional,
        Sequence,
        StrDict,
        Union,
    )

    _Molecules = Sequence[Molecule]


class FractionInStockScorer(Scorer):
    """Class for scoring nodes based on the fraction in stock"""

    scorer_name = "fraction in stock"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        # This is necessary because config should not be optional for this scorer
        self._config: Configuration = config

    def _score_node(self, node: MctsNode) -> float:
        num_in_stock = np.sum(node.state.in_stock_list)
        num_molecules = len(node.state.mols)
        return float(num_in_stock) / float(num_molecules)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        leaves = list(tree.leafs())
        num_in_stock = sum(mol in self._config.stock for mol in leaves)
        num_molecules = len(leaves)
        return float(num_in_stock) / float(num_molecules)


class FractionInSourceStockScorer(Scorer):
    """Class for scoring nodes based on the fraction in stock"""

    scorer_name = "fraction in source"

    def __init__(
        self,
        config: Configuration,
        source_stocks: List[str],
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        # This is necessary because config should not be optional for this scorer
        self._config: Configuration = config
        self._source_stocks = source_stocks
        self.scorer_name = f"fraction in {', '.join(source_stocks)}"

    def _calc_node_scores_for_source(self, node_mols: Sequence[Molecule]) -> float:
        num_in_stock = 0
        for mol in node_mols:
            availability = self._config.stock.availability_list(mol)
            sources = [
                source for source in availability if source in self._source_stocks
            ]
            if len(sources) > 0:
                num_in_stock += 1

        num_molecules = len(node_mols)
        return float(num_in_stock) / float(num_molecules)

    def _score_node(self, node: MctsNode) -> float:
        return self._calc_node_scores_for_source(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        leaves = list(tree.leafs())
        return self._calc_node_scores_for_source(leaves)


class FractionOfIntermediatesInStockScorer(Scorer):
    """Class for scoring nodes based on the fraction of intermediates in stock"""

    def __init__(
        self,
        config: Configuration,
        stock_name: str,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._stock: Stock = config.stock[stock_name]
        self._stock_name = stock_name
        self.scorer_name = f"fraction of intermediates in {self._stock_name}"

    def _fraction_in_stock(self, mols: Sequence[Molecule]) -> float:
        if len(mols) == 0:
            return 1.0
        return sum(mol in self._stock for mol in mols) / len(mols)

    def _score_node(self, node: MctsNode) -> int:
        reactions = node.actions_to()
        # Find intermediates by extacting the product from all reactions, except the
        # first reaction where the product is the root molecule
        intermediates = [reaction.mol for reaction in reactions[1::]]
        intermediates += node.state.expandable_mols
        return self._fraction_in_stock(intermediates)

    def _score_reaction_tree(self, tree: ReactionTree) -> int:
        tree_mols = list(tree.molecules())
        in_stock = get_node_attributes(tree.graph, "in_stock")
        mols = [mol for mol in tree_mols if not in_stock[mol] and mol is not tree.root]
        return self._fraction_in_stock(mols)


class NumberOfPrecursorsScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in a node or route"""

    scorer_name = "number of pre-cursors"

    def __init__(
        self,
        config: Optional[Configuration] = None,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._reverse_order = False

    def _score_node(self, node: MctsNode) -> float:
        return len(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.leafs()))


class NumberOfPrecursorsInStockScorer(Scorer):
    """Class for scoring nodes based on the number of pre-cursors in stock a
    node or route"""

    scorer_name = "number of pre-cursors in stock"

    def __init__(
        self, config: Configuration, scaler_params: Optional[StrDict] = None
    ) -> None:
        super().__init__(config, scaler_params)
        self._stock = config.stock

    def _score_node(self, node: MctsNode) -> float:
        return len([mol for mol in node.state.mols if mol in self._stock])

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len([mol for mol in tree.leafs() if mol in self._stock])


class PriceSumScorer(Scorer):
    """Scorer that sums the prices of all pre-cursors"""

    scorer_name = "sum of prices"

    def __init__(
        self,
        config: Configuration,
        scaler_params: Optional[StrDict] = None,
        default_cost: float = 1.0,
        not_in_stock_multiplier: int = 10,
    ) -> None:
        super().__init__(config, scaler_params)
        self._config: Configuration = config
        self.default_cost = default_cost
        self.not_in_stock_multiplier = not_in_stock_multiplier
        self._reverse_order = False

    def _calculate_leaf_costs(
        self, leafs: Union[Sequence[Molecule], Iterable[Molecule]]
    ) -> dict:
        costs = {}
        for mol in leafs:
            if mol not in self._config.stock:
                continue
            try:
                cost = self._config.stock.price(mol)
            except StockException:
                costs[mol] = self.default_cost
            else:
                costs[mol] = cost

        max_cost = max(costs.values()) if costs else self.default_cost
        return defaultdict(lambda: max_cost * self.not_in_stock_multiplier, costs)

    def _score_node(self, node: MctsNode) -> float:
        leaf_costs = self._calculate_leaf_costs(node.state.mols)
        return sum(leaf_costs[mol] for mol in node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        leaf_costs = self._calculate_leaf_costs(tree.leafs())
        return sum(leaf_costs[leaf] for leaf in tree.leafs())


class StockAvailabilityScorer(Scorer):
    """
    Scorer that computes score based on the stock availability of the starting material

    The score is calculated as a product of a stock score per starting material. The stock
    score for each molecule is based on the source of the stock, or a default value if the
    molecule is not in stock. The `other_source_score` parameter can be used to
    distinguish between "not in stock" and "not in the specificed sources" cases.
    """

    def __init__(
        self,
        config: Configuration,
        source_score: StrDict,
        default_score: float = 0.1,
        other_source_score: Optional[float] = None,
    ) -> None:
        super().__init__(config)
        assert self._config is not None
        self.source_score = source_score
        self.default_score = default_score
        self.other_source_score = other_source_score

    def __repr__(self) -> str:
        return "stock availability"

    def _calculate_leaf_costs(
        self, leafs: Union[Sequence[Molecule], Iterable[Molecule]]
    ) -> float:
        assert self._config is not None
        prod = 1.0
        for mol in leafs:
            availability = self._config.stock.availability_list(mol)
            scores = [
                self.source_score[source]
                for source in availability
                if source in self.source_score
            ]
            if scores:
                prod *= max(scores)
            elif self.other_source_score and mol in self._config.stock:
                prod *= self.other_source_score
            else:
                prod *= self.default_score
        return prod

    def _score_node(self, node: MctsNode) -> float:
        return self._calculate_leaf_costs(node.state.mols)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calculate_leaf_costs(tree.leafs())


class DeltaSyntheticComplexityScorer(Scorer):
    """
    Class for scoring nodes based on the delta-synthetic-complexity of the node
    and its parent 'horizon' steps up in the tree.

    :param config: the configuration the tree search
    :param sc_score_model: the path to the SCScore model
    :param scaler_params: the parameter settings of the scaler, defaults to max-min between -1.5 and 4
    :param horizon: the number of steps backwards to look for parent molecule
    """

    scorer_name: str = "delta-SC score"

    def __init__(
        self,
        config: Configuration,
        sc_score_model: str,
        scaler_params: Optional[StrDict] = None,
        horizon: int = 3,
    ) -> None:
        if scaler_params is None:
            scaler_params = {
                "name": "min_max",
                "min_val": -1.5,
                "max_val": 4,
                "reverse": False,
            }
        super().__init__(config, scaler_params)
        ## This is necessary because config should not be optional for this scorer
        self._config: Configuration = config

        self.horizon = horizon
        self._model = SCScore(sc_score_model)

    def sc_deltas(self, mols: _Molecules, parents: _Molecules) -> Sequence[float]:
        """
        Calculate delta SC-score among the list of mol-parent pairs.

        :params mols: the leaves of the tree
        :param parents: the parent of the leaves
        :returns: the pair-wise difference in SCScore
        """
        delta_sc_scores = []
        for mol, parent in zip(mols, parents):
            mol.sanitize()
            parent.sanitize()

            sc_score_mol = self._model(mol.rd_mol)
            sc_score_parent = self._model(parent.rd_mol)

            delta_sc_score = sc_score_parent - sc_score_mol
            delta_sc_scores.append(delta_sc_score)
        return delta_sc_scores

    def _get_parent_from_reaction_tree(
        self, tree: ReactionTree, mol: UniqueMolecule, horizon: int
    ) -> UniqueMolecule:
        horizon = min(horizon, tree.depth(mol) // 2)

        parent = mol
        # Step up in the tree, first a reaction, then a molecule
        for _ in range(horizon):
            parent = tree.parent_molecule(parent)

        return parent

    def _get_parent_from_tree_molecule(
        self, mol: TreeMolecule, horizon: int
    ) -> TreeMolecule:
        horizon = min(horizon, mol.transform)

        parent = mol
        # Step up in the tree via parent nodes
        for _ in range(horizon):
            grandparent = parent.parent
            if not grandparent:
                break
            parent = grandparent

        return parent

    def _score_node(self, node: MctsNode) -> float:
        expandable_mols = node.state.expandable_mols
        if not expandable_mols:
            expandable_mols = list(node.state.mols)

        parent_mols = [
            self._get_parent_from_tree_molecule(mol, self.horizon)
            for mol in expandable_mols
        ]
        delta_sc_scores = self.sc_deltas(expandable_mols, parent_mols)
        return min(delta_sc_scores)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        expandable_mols = [node for node in tree.leafs() if not tree.in_stock(node)]
        if not expandable_mols:
            expandable_mols = list(tree.leafs())

        parent_mols = [
            self._get_parent_from_reaction_tree(tree, node, self.horizon)
            for node in expandable_mols
        ]
        delta_sc_scores = self.sc_deltas(expandable_mols, parent_mols)
        return min(delta_sc_scores)
