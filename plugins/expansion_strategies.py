""" Module containing classes that implements different expansion policy strategies
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import requests
from requests.exceptions import ConnectionError  # pylint: disable=redefined-builtin

from aizynthfinder.chem import SmilesBasedRetroReaction
from aizynthfinder.context.policy import ExpansionStrategy
from aizynthfinder.utils.math import softmax

if TYPE_CHECKING:
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Dict, List, Sequence, Tuple

try:
    from ssbenchmark.model_zoo import ModelZoo
except ImportError:
    HAS_MODELZOO = False
else:
    HAS_MODELZOO = True


class ChemformerBasedExpansionStrategy(ExpansionStrategy):
    """
    A template-free expansion strategy that will return `SmilesBasedRetroReaction` objects upon expansion.
    It is based on calls to a REST API to the Chemformer model

    :param key: the key or label
    :param config: the configuration of the tree search
    :param url: the URL to the REST API
    :param ntrials: how many time to try a REST request
    """

    _required_kwargs = ["url"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        super().__init__(key, config, **kwargs)

        self._url: str = kwargs["url"]
        self._ntrials = kwargs.get("ntrials", 3)
        self._cache: Dict[str, Tuple[Sequence[str], Sequence[float]]] = {}

    # pylint: disable=R0914
    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Sequence[TreeMolecule] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules that are sent to
                                  the expansion model but for which predictions are not returned
        :return: the actions and the priors of those actions
        """
        possible_actions = []
        priors = []

        cache_molecules = cache_molecules or []
        self._update_cache(molecules + cache_molecules)

        for mol in molecules:
            try:
                output_smiles, probs = self._cache[mol.inchi_key]
            except KeyError:
                continue

            priors.extend(probs)
            for idx, reactants_str in enumerate(output_smiles):
                metadata = {}
                metadata["policy_probability"] = float(probs[idx])
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key
                possible_actions.append(
                    SmilesBasedRetroReaction(
                        mol, metadata=metadata, reactants_str=reactants_str
                    )
                )
        return possible_actions, priors  # type: ignore

    def reset_cache(self) -> None:
        """Reset the prediction cache"""
        self._cache = {}

    def _update_cache(self, molecules: Sequence[TreeMolecule]) -> None:
        pred_inchis = []
        smiles_list = []
        for molecule in molecules:
            if molecule.inchi_key in self._cache or molecule.inchi_key in pred_inchis:
                continue
            smiles_list.append(molecule.smiles)
            pred_inchis.append(molecule.inchi_key)

        if not pred_inchis:
            return

        for _ in range(self._ntrials):
            try:
                ret = requests.post(self._url, json=smiles_list)
            except ConnectionError:
                continue
            if ret.status_code == requests.codes.ok:
                break

        if ret.status_code != requests.codes.ok:
            self._logger.debug(
                f"Failed to retrieve results from Chemformer model: {ret.content}"
            )
            return

        predictions = ret.json()
        for prediction, inchi in zip(predictions, pred_inchis):
            self._cache[inchi] = (prediction["output"], softmax(prediction["lhs"]))


class ModelZooExpansionStrategy(ExpansionStrategy):
    """
    An expansion strategy that uses a single step model to operate on a Smiles-level
    of abstraction

    :param key: the key or label of the single step model
    :param config: the configuration of the tree search
    :param module_path: the path to the external model

    :raises ImportError: if ssbenchmark has not been installed.
    """

    _required_kwargs = ["module_path"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        if not HAS_MODELZOO:
            raise ImportError(
                "Cannot use this expansion strategy as it seems like "
                "ssbenchmark is not installed."
            )

        super().__init__(key, config, **kwargs)
        module_path = kwargs["module_path"]
        gpu_mode = kwargs.pop("use_gpu", False)
        model_params = dict(kwargs.pop("params", {}))

        self.model_zoo = ModelZoo(key, module_path, gpu_mode, model_params)
        self._cache: Dict[str, Tuple[Sequence[str], Sequence[float]]] = {}

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Sequence[TreeMolecule] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies
        and given cutoffs.

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules that are sent to the expansion
            model but for which predictions are not returned

        :return: the actions and the priors of those actions
        """
        possible_actions = []
        priors = []
        cache_molecules = cache_molecules or []
        self._update_cache(molecules + cache_molecules)

        for mol in molecules:
            output_smiles, probs = self._cache[mol.inchi_key]
            priors.extend(probs)

            for idx, move in enumerate(output_smiles):
                metadata = {}
                metadata["reaction"] = move
                metadata["policy_probability"] = float(probs[idx].round(4))
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key

                possible_actions.append(
                    SmilesBasedRetroReaction(mol, reactants_str=move, metadata=metadata)
                )

        return possible_actions, priors

    def _update_cache(self, molecules: Sequence[TreeMolecule]) -> None:
        pred_inchis = []
        smiles_list = []

        for molecule in molecules:
            if molecule.inchi_key in self._cache or molecule.inchi_key in pred_inchis:
                continue
            smiles_list.append(molecule.smiles)
            pred_inchis.append(molecule.inchi_key)

        if not pred_inchis:
            return

        pred_reactants, pred_priors = (
            np.array(item) for item in self.model_zoo.model_call(smiles_list)
        )

        for reactants, priors, inchi in zip(pred_reactants, pred_priors, pred_inchis):
            self._cache[inchi] = (reactants, priors)
