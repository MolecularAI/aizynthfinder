""" Module containing classes used to score the reaction routes.
These scores are only based on reactions in the route.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from rxnutils.routes.scoring import reaction_class_rank_score

from aizynthfinder.chem import TreeMolecule
from aizynthfinder.context.scoring.scorers_base import Scorer, make_rxnutils_route
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.search.mcts import MctsNode

if TYPE_CHECKING:
    from aizynthfinder.chem import FixedRetroReaction, RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import (
        List,
        Optional,
        Sequence,
        StrDict,
        Union,
        Tuple,
    )


class MaxTransformScorer(Scorer):
    """Class for scoring nodes based on the maximum transform"""

    scorer_name = "max transform"

    def _score_node(self, node: MctsNode) -> float:
        return node.state.max_transforms

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        mols = [
            TreeMolecule(
                parent=None, transform=tree.depth(leaf) // 2, smiles=leaf.smiles
            )
            for leaf in tree.leafs()
        ]
        return max(mol.transform for mol in mols)


class NumberOfReactionsScorer(Scorer):
    """Class for scoring nodes based on the number of reaction it took to get to a node"""

    scorer_name = "number of reactions"

    def __init__(
        self,
        config: Optional[Configuration] = None,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)
        self._reverse_order = False

    def _score_node(self, node: MctsNode) -> float:
        reactions = node.actions_to()
        return len(reactions)

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return len(list(tree.reactions()))


class AverageTemplateOccurrenceScorer(Scorer):
    """Class for scoring the nodes based on the average occurrence of the
    templates used to get to a node"""

    scorer_name = "average template occurrence"

    def _calc_average(
        self, reactions: Sequence[Union[FixedRetroReaction, RetroReaction]]
    ) -> float:
        if not reactions:
            return 0.0
        occurrences = [self._get_occurrence(reaction) for reaction in reactions]
        return sum(occurrences) / len(reactions)

    def _score_node(self, node: MctsNode) -> float:
        return self._calc_average(node.actions_to())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calc_average(list(tree.reactions()))

    @staticmethod
    def _get_occurrence(reaction: Union[FixedRetroReaction, RetroReaction]) -> int:
        return reaction.metadata.get(
            "library_occurrence", reaction.metadata.get("library_occurence", 0)
        )


class ReactionClassMembershipScorer(Scorer):
    """
    Scorer that checks if the reaction classes are in a specified set

    The score is calculated as product over each reaction. For each reaction
    the reaction classification is checked if it is in a given list of classes.
    """

    def __init__(
        self,
        config: Configuration,
        reaction_class_set: Union[str, List[str]],
        in_set_score: float = 1.0,
        not_in_set_score: float = 0.1,
    ) -> None:
        super().__init__(config)

        if isinstance(reaction_class_set, str):
            with open(reaction_class_set, "r") as file_id:
                self.reaction_class_set = [line.strip() for line in file_id.readlines()]
        else:
            self.reaction_class_set = reaction_class_set

        self.in_set_score = in_set_score
        self.not_in_set_score = not_in_set_score

    def __repr__(self) -> str:
        return "reaction class membership"

    def _calc_product(
        self, reactions: Sequence[Union[FixedRetroReaction, RetroReaction]]
    ) -> float:
        if not reactions:
            return 1.0
        prod = 1.0
        for reaction in reactions:
            membership = self._get_membership(reaction)
            prod *= self.in_set_score if membership else self.not_in_set_score
        return prod

    def _get_membership(
        self, reaction: Union[FixedRetroReaction, RetroReaction]
    ) -> bool:
        classification = reaction.metadata.get("classification")
        if not classification:
            return 0.0 in self.reaction_class_set
        return classification.split(" ")[0] in self.reaction_class_set

    def _score_node(self, node: MctsNode) -> float:
        return self._calc_product(node.actions_to())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        return self._calc_product(list(tree.reactions()))


class ReactionClassRankScorer(Scorer):
    """Class for scoring nodes and reaction trees by a weighted product of the ranks
    of the reaction classes in the route.

    Wrapper for `reaction_class_rank_score` from `rxnutils`

    :param config: the configuration the tree search
    :param class_ranks_path: the path to a file with reaction class ranks
    :param preferred_classes_path: the path to a file with preferred reaction classes
    :param non_preferred_factor: a penalty factor for non-preferred class
    :param scaler_params: the parameter settings of the scaler
    """

    scorer_name = "reaction class-rank score"

    def __init__(
        self,
        config: Configuration,
        class_ranks_path: str,
        preferred_classes_path: str,
        non_preferred_factor: float = 0.25,
        scaler_params: Optional[StrDict] = None,
    ) -> None:
        super().__init__(config, scaler_params)

        data = pd.read_csv(class_ranks_path, sep=",")
        self._reaction_class_ranks = dict(
            zip(data["reaction_class"], data["rank_score"])
        )

        with open(preferred_classes_path, "r") as fileobj:
            self._preferred_classes = fileobj.read().splitlines()

        self._non_preferred_factor = non_preferred_factor

    def _score_node(self, node: MctsNode) -> float:
        # All rxnutils scorers are based on reaction trees
        return self._score_reaction_tree(node.to_reaction_tree())

    def _score_reaction_tree(self, tree: ReactionTree) -> float:
        route = make_rxnutils_route(tree)
        return reaction_class_rank_score(
            route,
            self._reaction_class_ranks,
            self._preferred_classes,
            self._non_preferred_factor,
        )
