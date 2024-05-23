""" Module containing classes that implements different expansion policy strategies
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from aizynthfinder.chem import SmilesBasedRetroReaction, TemplatedRetroReaction
from aizynthfinder.context.policy.utils import _make_fingerprint
from aizynthfinder.utils.exceptions import PolicyException
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.models import load_model

if TYPE_CHECKING:
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        List,
        Optional,
        Sequence,
        StrDict,
        Tuple,
    )


class ExpansionStrategy(abc.ABC):
    """
    A base class for all expansion strategies.

    The strategy can be used by either calling the `get_actions` method
    of by calling the instantiated class with a list of molecule.

    .. code-block::

        expander = MyExpansionStrategy("dummy", config)
        actions, priors = expander.get_actions(molecules)
        actions, priors = expander(molecules)

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    _required_kwargs: List[str] = []

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise PolicyException(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self._config = config
        self._logger = logger()
        self.key = key

    def __call__(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        return self.get_actions(molecules, cache_molecules)

    @abc.abstractmethod
    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
                                  policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """

    def reset_cache(self) -> None:
        """Reset the prediction cache"""


class MultiExpansionStrategy(ExpansionStrategy):
    """
    A base class for combining multiple expansion strategies.

    The strategy can be used by either calling the `get_actions` method
    or by calling the instantiated class with a list of molecules.

    :ivar expansion_strategy_keys: the keys of the selected expansion strategies
    :ivar additive_expansion: a conditional setting to specify whether all the actions
        and priors of the selected expansion strategies should be combined or not.
        Defaults to False.
    :ivar expansion_strategy_weights: a list of weights for each expansion strategy.
        The weights should sum to one. Exception is the default, where unity weight
        is associated to each strategy.

    :param key: the key or label
    :param config: the configuration of the tree search
    :param expansion_strategies: the keys of the selected expansion strategies. All keys
        of the selected expansion strategies must exist in the expansion policies listed
        in config
    """

    _required_kwargs = ["expansion_strategies"]

    def __init__(
        self,
        key: str,
        config: Configuration,
        **kwargs: Any,
    ) -> None:
        super().__init__(key, config, **kwargs)
        self._config = config
        self._expansion_strategies: List[ExpansionStrategy] = []
        self.expansion_strategy_keys = kwargs["expansion_strategies"]

        self.cutoff_number = kwargs.get("cutoff_number")
        if self.cutoff_number:
            print(f"Setting multi-expansion cutoff_number: {self.cutoff_number}")

        self.expansion_strategy_weights = self._set_expansion_strategy_weights(kwargs)
        self.additive_expansion: bool = bool(kwargs.get("additive_expansion", False))
        self._logger.info(
            f"Multi-expansion strategy with policies: {self.expansion_strategy_keys}"
            f", and corresponding weights: {self.expansion_strategy_weights}"
        )

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies.

        The default implementation combines all the actions and priors of the
        selected expansion strategies into two lists respectively if the
        'additive_expansion' setting is set to True. This function can be overridden by
        a sub class to combine different expansion strategies in different ways.

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
            policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        :raises: PolicyException: if the policy isn't selected
        """
        expansion_strategies = self._get_expansion_strategies_from_config()

        all_possible_actions = []
        all_priors = []
        for expansion_strategy, expansion_strategy_weight in zip(
            expansion_strategies, self.expansion_strategy_weights
        ):
            possible_actions, priors = expansion_strategy.get_actions(
                molecules, cache_molecules
            )

            all_possible_actions.extend(possible_actions)
            if not self.additive_expansion and all_possible_actions:
                all_priors.extend(priors)
                break

            weighted_prior = [expansion_strategy_weight * p for p in priors]

            all_priors.extend(weighted_prior)

        all_possible_actions, all_priors = self._prune_actions(
            all_possible_actions, all_priors
        )
        return all_possible_actions, all_priors

    def _get_expansion_strategies_from_config(self) -> List[ExpansionStrategy]:
        if self._expansion_strategies:
            return self._expansion_strategies

        if not all(
            key in self._config.expansion_policy.items
            for key in self.expansion_strategy_keys
        ):
            raise ValueError(
                "The input expansion strategy keys must exist in the "
                "expansion policies listed in config"
            )
        self._expansion_strategies = [
            self._config.expansion_policy[key] for key in self.expansion_strategy_keys
        ]

        for expansion_strategy, weight in zip(
            self._expansion_strategies, self.expansion_strategy_weights
        ):
            if not getattr(expansion_strategy, "rescale_prior", True) and weight < 1:
                setattr(expansion_strategy, "rescale_prior", True)
                self._logger.info(
                    f"Enforcing {expansion_strategy.key}.rescale_prior=True"
                )
        return self._expansion_strategies

    def _prune_actions(
        self, actions: List[RetroReaction], priors: List[float]
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Prune the actions if a maximum number of actions is specified.

        :param actions: list of predicted actions
        :param priors: list of prediction probabilities
        :return: the top 'self.cutoff_number' actions and corresponding priors.
        """
        if not self.cutoff_number:
            return actions, priors

        sortidx = np.argsort(np.array(priors))[::-1].astype(int)
        priors = [priors[idx] for idx in sortidx[0 : self.cutoff_number]]
        actions = [actions[idx] for idx in sortidx[0 : self.cutoff_number]]
        return actions, priors

    def _set_expansion_strategy_weights(self, kwargs: StrDict) -> List[float]:
        """
        Set the weights of each expansion strategy using the input kwargs from config.
        The weights in the config should sum to one.
        If not set in the config file, the weights default to one for each strategy
        (for backwards compatibility).

        :param kwargs: input arguments to the MultiExpansionStrategy
        :raises: ValueError if weights from the config file do not sum to one.
        :return: a list of expansion strategy weights
        """
        if not "expansion_strategy_weights" in kwargs:
            return [1.0 for _ in self.expansion_strategy_keys]

        expansion_strategy_weights = kwargs["expansion_strategy_weights"]
        sum_weights = sum(expansion_strategy_weights)

        if sum_weights != 1:
            raise ValueError(
                "The expansion strategy weights in MultiExpansion should "
                "sum to one. -> "
                f"sum({expansion_strategy_weights})={sum_weights}."
            )

        return expansion_strategy_weights


class TemplateBasedExpansionStrategy(ExpansionStrategy):
    """
    A template-based expansion strategy that will return `TemplatedRetroReaction` objects upon expansion.

    :ivar template_column: the column in the template file that contains the templates
    :ivar cutoff_cumulative: the accumulative probability of the suggested templates
    :ivar cutoff_number: the maximum number of templates to returned
    :ivar use_rdchiral: a boolean to apply templates with RDChiral
    :ivar use_remote_models: a boolean to connect to remote TensorFlow servers
    :ivar rescale_prior: a boolean to apply softmax to the priors
    :ivar chiral_fingerprints: if True will base expansion on chiral fingerprint
    :ivar mask: a boolean vector of masks for the reaction templates. The length of the vector should be equal to the
        number of templates. It is set to None if no mask file is provided as input.

    :param key: the key or label
    :param config: the configuration of the tree search
    :param model: the source of the policy model
    :param template: the path to a HDF5 file with the templates
    :raises PolicyException: if the length of the model output vector is not same as the
        number of templates
    """

    _required_kwargs = [
        "model",
        "template",
    ]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        super().__init__(key, config, **kwargs)

        source = kwargs["model"]
        templatefile = kwargs["template"]
        maskfile: str = kwargs.get("mask", "")
        self.template_column: str = kwargs.get("template_column", "retro_template")
        self.cutoff_cumulative: float = float(kwargs.get("cutoff_cumulative", 0.995))
        self.cutoff_number: int = int(kwargs.get("cutoff_number", 50))
        self.use_rdchiral: bool = bool(kwargs.get("use_rdchiral", True))
        self.use_remote_models: bool = bool(kwargs.get("use_remote_models", False))
        self.rescale_prior: bool = bool(kwargs.get("rescale_prior", False))
        self.chiral_fingerprints = bool(kwargs.get("chiral_fingerprints", False))

        self._logger.info(
            f"Loading template-based expansion policy model from {source} to {self.key}"
        )
        self.model = load_model(source, self.key, self.use_remote_models)

        self._logger.info(f"Loading templates from {templatefile} to {self.key}")
        if templatefile.endswith(".csv.gz") or templatefile.endswith(".csv"):
            self.templates: pd.DataFrame = pd.read_csv(
                templatefile, index_col=0, sep="\t"
            )
        else:
            self.templates = pd.read_hdf(templatefile, "table")

        self.mask: Optional[np.ndarray] = (
            self._load_mask_file(maskfile) if maskfile else None
        )

        if hasattr(self.model, "output_size") and len(self.templates) != self.model.output_size:  # type: ignore
            raise PolicyException(
                f"The number of templates ({len(self.templates)}) does not agree with the "  # type: ignore
                f"output dimensions of the model ({self.model.output_size})"
            )
        self._cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
                                  policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """

        possible_actions = []
        priors: List[float] = []
        cache_molecules = cache_molecules or []
        self._update_cache(list(molecules) + list(cache_molecules))

        for mol in molecules:
            probable_transforms_idx, probs = self._cache[mol.inchi_key]
            possible_moves = self.templates.iloc[probable_transforms_idx]
            if self.rescale_prior:
                probs /= probs.sum()
            priors.extend(probs)
            for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
                metadata = dict(move)
                del metadata[self.template_column]
                metadata["policy_probability"] = float(probs[idx].round(4))
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key
                metadata["template_code"] = move_index
                metadata["template"] = move[self.template_column]
                possible_actions.append(
                    TemplatedRetroReaction(
                        mol,
                        smarts=move[self.template_column],
                        metadata=metadata,
                        use_rdchiral=self.use_rdchiral,
                    )
                )
        return possible_actions, priors  # type: ignore

    def reset_cache(self) -> None:
        """Reset the prediction cache"""
        self._cache = {}

    def _cutoff_predictions(self, predictions: np.ndarray) -> np.ndarray:
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        if self.mask is not None:
            predictions[~self.mask] = 0
        sortidx = np.argsort(predictions)[::-1]
        cumsum: np.ndarray = np.cumsum(predictions[sortidx])
        if any(cumsum >= self.cutoff_cumulative):
            maxidx = int(np.argmin(cumsum < self.cutoff_cumulative))
        else:
            maxidx = len(cumsum)
        maxidx = min(maxidx, self.cutoff_number) or 1
        return sortidx[:maxidx]

    def _load_mask_file(self, maskfile: str) -> np.ndarray:
        self._logger.info(f"Loading masking of templates from {maskfile} to {self.key}")
        mask = np.load(maskfile)["arr_0"]
        if len(mask) != len(self.templates):
            raise PolicyException(
                f"The number of masks {len(mask)} does not match the number of templates {len(self.templates)}"
            )
        return mask

    def _update_cache(self, molecules: Sequence[TreeMolecule]) -> None:
        pred_inchis = []
        fp_list = []
        for molecule in molecules:
            if molecule.inchi_key in self._cache or molecule.inchi_key in pred_inchis:
                continue
            fp_list.append(
                _make_fingerprint(molecule, self.model, self.chiral_fingerprints)
            )
            pred_inchis.append(molecule.inchi_key)

        if not pred_inchis:
            return

        pred_list = np.asarray(self.model.predict(np.vstack(fp_list)))
        for pred, inchi in zip(pred_list, pred_inchis):
            probable_transforms_idx = self._cutoff_predictions(pred)
            self._cache[inchi] = (
                probable_transforms_idx,
                pred[probable_transforms_idx],
            )


class TemplateBasedDirectExpansionStrategy(TemplateBasedExpansionStrategy):
    """
    A template-based expansion strategy that will return `SmilesBasedRetroReaction` objects upon expansion
    by directly applying the template

    :param key: the key or label
    :param config: the configuration of the tree search
    :param source: the source of the policy model
    :param templatefile: the path to a HDF5 file with the templates
    :raises PolicyException: if the length of the model output vector is not same as the number of templates
    """

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
            policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """
        possible_actions = []
        priors = []

        super_actions, super_priors = super().get_actions(molecules, cache_molecules)
        for templated_action, prior in zip(super_actions, super_priors):
            for reactants in templated_action.reactants:
                reactants_str = ".".join(mol.smiles for mol in reactants)
                new_action = SmilesBasedRetroReaction(
                    templated_action.mol,
                    metadata=templated_action.metadata,
                    reactants_str=reactants_str,
                )
                possible_actions.append(new_action)
                priors.append(prior)

        return possible_actions, priors  # type: ignore
