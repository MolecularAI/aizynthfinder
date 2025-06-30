""" Module containing classes that implements different filter policy strategies
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

import numpy as np
from rdkit import Chem

from aizynthfinder.chem import TemplatedRetroReaction
from aizynthfinder.context.policy.utils import _make_fingerprint
from aizynthfinder.utils.bonds import BrokenBonds
from aizynthfinder.utils.exceptions import PolicyException, RejectionException
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.models import load_model

if TYPE_CHECKING:
    from aizynthfinder.chem.reaction import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Any, Dict, List, Optional, Tuple


class FilterStrategy(abc.ABC):
    """
    A base class for all filter strategies.

    The filter can be applied by either calling the `apply` method
    of by calling the instantiated class with a reaction.

    .. code-block::

        filter = MyFilterStrategy("dummy", config)
        filter.apply(reaction)
        filter(reaction)

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    _required_kwargs: List[str] = []

    def __init__(self, key: str, config: Configuration, **kwargs: Any) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise PolicyException(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self._config = config
        self._logger = logger()
        self.key = key

    def __call__(self, reaction: RetroReaction) -> None:
        self.apply(reaction)

    @abc.abstractmethod
    def apply(self, reaction: RetroReaction) -> None:
        """
        Apply the filter on the reaction. If the reaction
        should be rejected a `RejectionException` is raised

        :param reaction: the reaction to filter
        :raises: if the reaction should be rejected.
        """


class BondFilter(FilterStrategy):
    """
    Check if focussed bonds to freeze stay frozen in a reaction.

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    def __init__(self, key: str, config: Configuration, **kwargs: Any) -> None:
        super().__init__(key, config, **kwargs)

        self._freeze_bonds = config.search.freeze_bonds
        self._broken_bonds = BrokenBonds(self._freeze_bonds)
        self._logger.info(
            f"Loading bond filter to {key} with {len(self._freeze_bonds)} "
            "bonds to freeze"
        )

    def apply(self, reaction: RetroReaction) -> None:
        broken_frozen_bonds = self._broken_bonds(reaction)
        if len(broken_frozen_bonds) > 0:
            raise RejectionException(
                f"{reaction} was filtered out as the focussed bonds "
                f"'{broken_frozen_bonds}' were found to be broken in the reaction"
            )


class QuickKerasFilter(FilterStrategy):
    """
    Filter quick-filter trained on artificial negative data

    :ivar use_remote_models: a boolean to connect to remote TensorFlow servers. Defaults
        to False.
    :ivar filter_cutoff: the cut-off value

    :param key: the key or label
    :param config: the configuration of the tree search
    :param model: the source of the policy model
    """

    _required_kwargs: List[str] = ["model"]

    def __init__(self, key: str, config: Configuration, **kwargs: Any) -> None:
        super().__init__(key, config, **kwargs)
        source = kwargs["model"]
        # self.settings = self._config.filter_settings
        self._logger.info(f"Loading filter policy model from {source} to {key}")
        self.use_remote_models: bool = bool(kwargs.get("use_remote_models", False))
        self.model = load_model(source, key, self.use_remote_models)
        self._prod_fp_name = kwargs.get("prod_fp_name", "input_1")
        self._rxn_fp_name = kwargs.get("rxn_fp_name", "input_2")
        self._exclude_from_policy: List[str] = kwargs.get("exclude_from_policy", [])
        self.filter_cutoff: float = float(kwargs.get("filter_cutoff", 0.05))

    def apply(self, reaction: RetroReaction) -> None:
        if reaction.metadata.get("policy_name", "") in self._exclude_from_policy:
            return

        feasible, prob = self.feasibility(reaction)
        if not feasible:
            raise RejectionException(f"{reaction} was filtered out with prob {prob}")

    def feasibility(self, reaction: RetroReaction) -> Tuple[bool, float]:
        """
        Computes if a given reaction is feasible by given
        the reaction fingerprint to a network model

        :param reaction: the reaction to query
        :return: if the reaction is feasible
        """
        if not reaction.reactants:
            return False, 0.0

        prob = self._predict(reaction)
        feasible = prob >= self.filter_cutoff
        return feasible, prob

    def _predict(self, reaction: RetroReaction) -> float:
        prod_fp, rxn_fp = self._reaction_to_fingerprint(reaction, self.model)
        kwargs = {self._prod_fp_name: prod_fp, self._rxn_fp_name: rxn_fp}
        return self.model.predict(prod_fp, rxn_fp, **kwargs)[0][0]

    @staticmethod
    def _reaction_to_fingerprint(
        reaction: RetroReaction, model: Any
    ) -> Tuple[np.ndarray, np.ndarray]:
        rxn_fp = _make_fingerprint(reaction, model)
        prod_fp = _make_fingerprint(reaction.mol, model)
        return prod_fp, rxn_fp


class ReactantsCountFilter(FilterStrategy):
    """
    Check that the number of reactants is was expected from the template

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    def __init__(self, key: str, config: Configuration, **kwargs: Any) -> None:
        super().__init__(key, config, **kwargs)
        self._logger.info(f"Loading reactants count filter to {key}")

    def apply(self, reaction: RetroReaction) -> None:
        if not isinstance(reaction, TemplatedRetroReaction):
            raise ValueError(
                "Reactants count filter can only be used on templated retro reaction "
            )

        reactants = reaction.reactants[reaction.index]
        if len(reactants) > reaction.rd_reaction.GetNumProductTemplates():
            raise RejectionException(
                f"{reaction} was filtered out because number of reactants disagree with the template"
            )


class FrozenSubstructureFilter(FilterStrategy):
    """
    Filter for rejecting reactions that break up substructures

    :param key: the key or label
    :param config: the configuration of the tree search
    :param smarts_list: the SMARTS patterns of the sub-structures
    """

    _required_kwargs: List[str] = ["smarts_list"]

    def __init__(self, key: str, config: Configuration, **kwargs: Any) -> None:
        super().__init__(key, config, **kwargs)
        self._smarts_list = kwargs.get("smarts_list", [])
        self._mol_lists = [Chem.MolFromSmarts(smarts) for smarts in self._smarts_list]
        self._logger.info(
            f"Loading frozen substructure filter to {key} with {len(self._mol_lists)} substructures"
        )

    def apply(self, reaction: RetroReaction) -> None:
        for mol in self._mol_lists:
            # If it did not exists in the product, we cannot expect it to be present in the reactants
            if not reaction.mol.rd_mol.HasSubstructMatch(mol):
                continue

            found = False
            for reactant in reaction.reactants[reaction.index]:
                if reactant.rd_mol.HasSubstructMatch(mol):
                    found = True
                    break

            if not found:
                raise RejectionException(
                    f"{reaction} was filtered out because of broken substructure: {Chem.MolToSmarts(mol)}"
                )


FILTER_STRATEGY_ALIAS = {
    "feasibility": "QuickKerasFilter",
    "quick_keras_filter": "QuickKerasFilter",
    "reactants_count": "ReactantsCountFilter",
    "frozen_substructure": "FrozenSubstructureFilter",
}
