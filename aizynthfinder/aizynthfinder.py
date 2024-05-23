""" Module containing a class that is the main interface the retrosynthesis tool.
"""
from __future__ import annotations

import time
from collections import defaultdict
from typing import TYPE_CHECKING

from tqdm import tqdm

from aizynthfinder.analysis import (
    RouteCollection,
    RouteSelectionArguments,
    TreeAnalysis,
)
from aizynthfinder.chem import FixedRetroReaction, Molecule, TreeMolecule
from aizynthfinder.context.config import Configuration
from aizynthfinder.context.policy import BondFilter
from aizynthfinder.context.scoring import BrokenBondsScorer, CombinedScorer
from aizynthfinder.reactiontree import ReactionTreeFromExpansion
from aizynthfinder.search.andor_trees import AndOrSearchTreeBase
from aizynthfinder.search.mcts import MctsSearchTree
from aizynthfinder.utils.exceptions import MoleculeException
from aizynthfinder.utils.loading import load_dynamic_class

# This must be imported first to setup logging for rdkit, tensorflow etc
from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.chem import RetroReaction
    from aizynthfinder.utils.type_utils import (
        Callable,
        Dict,
        List,
        Optional,
        StrDict,
        Tuple,
        Union,
    )


class AiZynthFinder:
    """
    Public API to the aizynthfinder tool

    If instantiated with the path to a yaml file or dictionary of settings
    the stocks and policy networks are loaded directly.
    Otherwise, the user is responsible for loading them prior to
    executing the tree search.

    :ivar config: the configuration of the search
    :ivar expansion_policy: the expansion policy model
    :ivar filter_policy: the filter policy model
    :ivar stock: the stock
    :ivar scorers: the loaded scores
    :ivar tree: the search tree
    :ivar analysis: the tree analysis
    :ivar routes: the top-ranked routes
    :ivar search_stats: statistics of the latest search

    :param configfile: the path to yaml file with configuration (has priority over configdict), defaults to None
    :param configdict: the config as a dictionary source, defaults to None
    """

    def __init__(
        self, configfile: Optional[str] = None, configdict: Optional[StrDict] = None
    ) -> None:
        self._logger = logger()

        if configfile:
            self.config = Configuration.from_file(configfile)
        elif configdict:
            self.config = Configuration.from_dict(configdict)
        else:
            self.config = Configuration()

        self.expansion_policy = self.config.expansion_policy
        self.filter_policy = self.config.filter_policy
        self.stock = self.config.stock
        self.scorers = self.config.scorers
        self.tree: Optional[Union[MctsSearchTree, AndOrSearchTreeBase]] = None
        self._target_mol: Optional[Molecule] = None
        self.search_stats: StrDict = dict()
        self.routes = RouteCollection([])
        self.analysis: Optional[TreeAnalysis] = None
        self._num_objectives = len(
            self.config.search.algorithm_config.get("search_rewards", [])
        )

    @property
    def target_smiles(self) -> str:
        """The SMILES representation of the molecule to predict routes on."""
        if not self._target_mol:
            return ""
        return self._target_mol.smiles

    @target_smiles.setter
    def target_smiles(self, smiles: str) -> None:
        self.target_mol = Molecule(smiles=smiles)

    @property
    def target_mol(self) -> Optional[Molecule]:
        """The molecule to predict routes on"""
        return self._target_mol

    @target_mol.setter
    def target_mol(self, mol: Molecule) -> None:
        self.tree = None
        self._target_mol = mol

    def build_routes(
        self,
        selection: Optional[RouteSelectionArguments] = None,
        scorer: Optional[Union[str, List[str]]] = None,
    ) -> None:
        """
        Build reaction routes

        This is necessary to call after the tree search has completed in order
        to extract results from the tree search.

        :param selection: the selection criteria for the routes
        :param scorer: a reference to the object used to score the nodes, can be a list
        :raises ValueError: if the search tree not initialized
        """
        self.analysis = self._setup_analysis(scorer=scorer)
        config_selection = RouteSelectionArguments(
            nmin=self.config.post_processing.min_routes,
            nmax=self.config.post_processing.max_routes,
            return_all=self.config.post_processing.all_routes,
        )
        self.routes = RouteCollection.from_analysis(
            self.analysis, selection or config_selection
        )

    def extract_statistics(self) -> StrDict:
        """Extracts tree statistics as a dictionary"""
        if not self.analysis:
            return {}
        stats = {
            "target": self.target_smiles,
            "search_time": self.search_stats["time"],
            "first_solution_time": self.search_stats.get("first_solution_time", 0),
            "first_solution_iteration": self.search_stats.get(
                "first_solution_iteration", 0
            ),
        }
        stats.update(self.analysis.tree_statistics())
        return stats

    def prepare_tree(self) -> None:
        """
        Setup the tree for searching

        :raises ValueError: if the target molecule was not set
        """
        if not self.target_mol:
            raise ValueError("No target molecule set")

        try:
            self.target_mol.sanitize()
        except MoleculeException:
            raise ValueError("Target molecule unsanitizable")

        self.stock.reset_exclusion_list()
        if (
            self.config.search.exclude_target_from_stock
            and self.target_mol in self.stock
        ):
            self.stock.exclude(self.target_mol)
            self._logger.debug("Excluding the target compound from the stock")

        if self.config.search.break_bonds or self.config.search.freeze_bonds:
            self._setup_focussed_bonds(self.target_mol)

        self._setup_search_tree()
        self.analysis = None
        self.routes = RouteCollection([])
        self.filter_policy.reset_cache()
        self.expansion_policy.reset_cache()

    def stock_info(self) -> StrDict:
        """
        Return the stock availability for all leaf nodes in all collected reaction trees

        The key of the return dictionary will be the SMILES string of the leaves,
        and the value will be the stock availability

        :return: the collected stock information.
        """
        if not self.analysis:
            return {}
        _stock_info = {}
        for tree in self.routes.reaction_trees:
            for leaf in tree.leafs():
                if leaf.smiles not in _stock_info:
                    _stock_info[leaf.smiles] = self.stock.availability_list(leaf)
        return _stock_info

    def tree_search(self, show_progress: bool = False) -> float:
        """
        Perform the actual tree search

        :param show_progress: if True, shows a progress bar
        :return: the time past in seconds
        """
        if not self.tree:
            self.prepare_tree()
        # This is for type checking, prepare_tree is creating it.
        assert self.tree is not None
        self.search_stats = {"returned_first": False, "iterations": 0}

        time0 = time.time()
        i = 1
        self._logger.debug("Starting search")
        time_past = time.time() - time0

        if show_progress:
            pbar = tqdm(total=self.config.search.iteration_limit, leave=False)

        while (
            time_past < self.config.search.time_limit
            and i <= self.config.search.iteration_limit
        ):
            if show_progress:
                pbar.update(1)
            self.search_stats["iterations"] += 1

            try:
                is_solved = self.tree.one_iteration()
            except StopIteration:
                break

            if is_solved and "first_solution_time" not in self.search_stats:
                self.search_stats["first_solution_time"] = time.time() - time0
                self.search_stats["first_solution_iteration"] = i

            if self.config.search.return_first and is_solved:
                self._logger.debug("Found first solved route")
                self.search_stats["returned_first"] = True
                break
            i = i + 1
            time_past = time.time() - time0

        if show_progress:
            pbar.close()
        time_past = time.time() - time0
        self._logger.debug("Search completed")
        self.search_stats["time"] = time_past
        return time_past

    def _setup_focussed_bonds(self, target_mol: Molecule) -> None:
        """
        Setup multi-objective scoring function with 'broken bonds'-scorer and
        add 'frozen bonds'-filter to filter policy.

        :param target_mol: the target molecule.
        """
        target_mol = TreeMolecule(smiles=target_mol.smiles, parent=None)

        bond_filter_key = "__finder_bond_filter"
        if self.config.search.freeze_bonds:
            if not target_mol.has_all_focussed_bonds(self.config.search.freeze_bonds):
                raise ValueError("Bonds in 'freeze_bond' must exist in target molecule")
            bond_filter = BondFilter(bond_filter_key, self.config)
            self.filter_policy.load(bond_filter)
            self.filter_policy.select(bond_filter_key, append=True)
        elif (
            self.filter_policy.selection
            and bond_filter_key in self.filter_policy.selection
        ):
            self.filter_policy.deselect(bond_filter_key)

        search_rewards = self.config.search.algorithm_config.get("search_rewards")
        if not search_rewards:
            return

        if self.config.search.break_bonds and "broken bonds" in search_rewards:
            if not target_mol.has_all_focussed_bonds(self.config.search.break_bonds):
                raise ValueError("Bonds in 'break_bonds' must exist in target molecule")
            self.scorers.load(BrokenBondsScorer(self.config))
            self._num_objectives = len(search_rewards)

    def _setup_search_tree(self) -> None:
        self._logger.debug(f"Defining tree root:  {self.target_smiles}")
        if self.config.search.algorithm.lower() == "mcts":
            self.tree = MctsSearchTree(
                root_smiles=self.target_smiles, config=self.config
            )
        else:
            cls = load_dynamic_class(self.config.search.algorithm)
            self.tree = cls(root_smiles=self.target_smiles, config=self.config)

    def _setup_analysis(
        self,
        scorer: Optional[Union[str, List[str]]],
    ) -> TreeAnalysis:
        """Configure TreeAnalysis

        :param scorer: a reference to the object used to score the nodes, can be a list
        :returns: the configured TreeAnalysis
        :raises ValueError: if the search tree not initialized
        """
        if not self.tree:
            raise ValueError("Search tree not initialized")

        if scorer is None:
            scorer_names = self.config.post_processing.route_scorers
            # If not defined, use the same scorer as the search rewards
            if not scorer_names:
                search_rewards = self.config.search.algorithm_config.get(
                    "search_rewards"
                )
                scorer_names = search_rewards if search_rewards else ["state score"]

        elif isinstance(scorer, str):
            scorer_names = [scorer]
        else:
            scorer_names = list(scorer)

        if "broken bonds" in scorer_names:
            # Add broken bonds scorer if required
            self.scorers.load(BrokenBondsScorer(self.config))

        scorers = [self.scorers[name] for name in scorer_names]

        if self.config.post_processing.scorer_weights:
            scorers = [
                CombinedScorer(
                    self.config,
                    scorer_names,
                    self.config.post_processing.scorer_weights,
                )
            ]

        return TreeAnalysis(self.tree, scorers)


class AiZynthExpander:
    """
    Public API to the AiZynthFinder expansion and filter policies

    If instantiated with the path to a yaml file or dictionary of settings
    the stocks and policy networks are loaded directly.
    Otherwise, the user is responsible for loading them prior to
    executing the tree search.

    :ivar config: the configuration of the search
    :ivar expansion_policy: the expansion policy model
    :ivar filter_policy: the filter policy model

    :param configfile: the path to yaml file with configuration (has priority over configdict), defaults to None
    :param configdict: the config as a dictionary source, defaults to None
    """

    def __init__(
        self, configfile: Optional[str] = None, configdict: Optional[StrDict] = None
    ) -> None:
        self._logger = logger()

        if configfile:
            self.config = Configuration.from_file(configfile)
        elif configdict:
            self.config = Configuration.from_dict(configdict)
        else:
            self.config = Configuration()

        self.expansion_policy = self.config.expansion_policy
        self.filter_policy = self.config.filter_policy
        self.stats: StrDict = {}

    def do_expansion(
        self,
        smiles: str,
        return_n: int = 5,
        filter_func: Optional[Callable[[RetroReaction], bool]] = None,
    ) -> List[Tuple[FixedRetroReaction, ...]]:
        """
        Do the expansion of the given molecule returning a list of
        reaction tuples. Each tuple in the list contains reactions
        producing the same reactants. Hence, nested structure of the
        return value is way to group reactions.

        If filter policy is setup, the probability of the reactions are
        added as metadata to the reaction.

        The additional filter functions makes it possible to do customized
        filtering. The callable should take as only argument a `RetroReaction`
        object and return True if the reaction can be kept or False if it should
        be removed.

        :param smiles: the SMILES string of the target molecule
        :param return_n: the length of the return list
        :param filter_func: an additional filter function
        :return: the grouped reactions
        """
        self.stats = {"non-applicable": 0}

        mol = TreeMolecule(parent=None, smiles=smiles)
        actions, _ = self.expansion_policy.get_actions([mol])
        results: Dict[Tuple[str, ...], List[FixedRetroReaction]] = defaultdict(list)
        for action in actions:
            reactants = action.reactants
            if not reactants:
                self.stats["non-applicable"] += 1
                continue
            if filter_func and not filter_func(action):
                continue
            for name in self.filter_policy.selection or []:
                if hasattr(self.filter_policy[name], "feasibility"):
                    _, feasibility_prob = self.filter_policy[name].feasibility(action)
                    action.metadata["feasibility"] = float(feasibility_prob)
                    break
            action.metadata["expansion_rank"] = len(results) + 1
            unique_key = tuple(sorted(mol.inchi_key for mol in reactants[0]))
            if unique_key not in results and len(results) >= return_n:
                continue
            rxn = next(ReactionTreeFromExpansion(action).tree.reactions())  # type: ignore
            results[unique_key].append(rxn)
        return [tuple(reactions) for reactions in results.values()]
