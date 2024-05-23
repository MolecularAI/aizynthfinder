""" Module containing classes that implements different expansion policy strategies
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import requests
from requests.exceptions import ConnectionError  # pylint: disable=redefined-builtin

from aizynthfinder.chem import SmilesBasedRetroReaction
from aizynthfinder.context.policy import ExpansionStrategy
from aizynthfinder.utils.bonds import BrokenBonds
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.math import softmax

if TYPE_CHECKING:
    from typing import Any, Optional
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import Dict, List, Sequence, Tuple, Union

try:
    from ssbenchmark.model_zoo import ModelZoo
except ImportError:
    HAS_MODELZOO = False
else:
    HAS_MODELZOO = True


class ChemformerBasedExpansionStrategy(ExpansionStrategy):
    """
    A template-free expansion strategy that will return `SmilesBasedRetroReaction`
    objects upon expansion.
    It is based on calls to a REST API to the Chemformer model

    :param key: the key or label
    :param config: the configuration of the tree search
    :param url: the URL to the REST API
    :param ntrials: how many times to try a REST request
    """

    _required_kwargs = ["url"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:

        super().__init__(key, config, **kwargs)

        self._ntrials = kwargs.get("ntrials", 3)
        self._n_beams = kwargs.get("n_beams", 10)
        self._model_url: str = kwargs["url"]

        self._cache: Dict[str, Tuple[Sequence[str], Sequence[float]]] = {}
        self._logger = logger()
        self._logger.info(
            f"Loading Chemformer with num predictions: {self._n_beams} to {self.key}"
        )

    # pylint: disable=R0914
    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Sequence[TreeMolecule] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions (list of SmilesBasedRetroReaction:s) of a set of molecules

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
            cache_key = self._cache_key(mol)
            try:
                predictions, probabilities = self._cache[cache_key]
            except KeyError:
                continue

            priors.extend(probabilities)
            for idx, reactants_smiles in enumerate(predictions):
                metadata = {}
                metadata["policy_probability"] = float(probabilities[idx])
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key
                possible_actions.append(
                    SmilesBasedRetroReaction(
                        mol, metadata=metadata, reactants_str=reactants_smiles
                    )
                )
        return possible_actions, priors  # type: ignore

    def reset_cache(self) -> None:
        """Reset the prediction cache"""
        self._cache = {}

    def _cache_key(self, molecule: TreeMolecule):
        return molecule.smiles

    def _get_predictions(
        self, model_input: Optional[Union[Dict[str, Any], List[str]]]
    ) -> Optional[requests.models.Response]:
        """Get predictions from model api."""
        if not model_input:
            return

        response = None
        for _ in range(self._ntrials):
            try:
                response = requests.post(
                    self._model_url,
                    json=model_input,
                    params={"n_beams": self._n_beams},
                )
            except ConnectionError:
                continue
            else:
                if response.status_code == requests.codes.ok:
                    break

        if not response:
            self._logger.debug(
                f"Failed to retrieve results from Chemformer model with url: "
                f"{self._model_url}, and input: {model_input}"
            )
            return

        if response.status_code != requests.codes.ok:
            self._logger.debug(
                f"Failed to retrieve results from Chemformer model: {response.content}"
            )
            return
        return response

    def _make_cache_item(self, prediction: Dict[str, Any]) -> Tuple[Any, ...]:
        return (prediction["output"], softmax(prediction["lhs"]))

    def _make_model_input(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[str], Union[Dict[str, Any], List[str]]]:
        """
        Construct input for the standard Chemformer model.
        :param molecules: a list of molecules
        :return: model input as a dictionary or None if no molecules to predict on
        """
        product_cache_keys = []
        input_smiles = []
        for molecule in molecules:

            cache_key = self._cache_key(molecule)

            if cache_key in self._cache or cache_key in product_cache_keys:
                continue

            input_smiles.append(molecule.smiles)
            product_cache_keys.append(cache_key)

        if not product_cache_keys:
            return [], []
        return product_cache_keys, input_smiles

    def _update_cache(self, molecules: Sequence[TreeMolecule], **kwargs) -> None:
        """
        Run retrosynthesis prediction on molecules which are not cached and update
        cache with the prediction outcomes.
        :param molecules: list of molecules
        :param bonds_to_break: for each molecule, list of bonds to disconnect.
        """
        cache_keys, model_input = self._make_model_input(molecules, **kwargs)
        response = self._get_predictions(model_input)

        if not response:
            return

        predictions = response.json()
        for prediction, cache_key in zip(predictions, cache_keys):
            self._cache[cache_key] = self._make_cache_item(prediction)


class DisconnectionAwareExpansionStrategy(ChemformerBasedExpansionStrategy):
    """
    A template-free expansion strategy that will return `SmilesBasedRetroReaction`
    objects upon expansion.
    It is based on calls to a REST API to the disconnection-aware Chemformer model.
    The disconnection-aware Chemformer is used to disconnect specific bonds when the user
    inputs bond constraints in the config file (config.search.break_bonds).
    :param key: the key or label
    :param config: the configuration of the tree search
    :param url: the URL to the REST API
    :param ntrials: how many times to try a REST request
    """

    _required_kwargs = ["url"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:

        super().__init__(key, config, **kwargs)

        if not config.search.break_bonds:
            raise ValueError(
                "Disconnection-aware Chemformer requires config.search.break_bonds to be defined."
            )

        self._precache_predictions = kwargs.get("precache_predictions", False)
        self._bonds_to_break = config.search.break_bonds
        self._disconnected_bonds = BrokenBonds(self._bonds_to_break)

        self._logger.info(
            f"Bonds to break constraints in Chemformer: {self._bonds_to_break}"
        )

    # pylint: disable=R0914
    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Sequence[TreeMolecule] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions (list of SmilesBasedRetroReaction:s) of a set of molecules
        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules that are sent to
                        the expansion model but for which predictions are not returned
        :return: the actions and the priors of those actions
        """

        possible_actions = []
        priors = []

        cache_molecules = cache_molecules or []

        bonds_to_break = []
        cache_bonds_to_break = []

        for mol in molecules:
            bonds_to_break.append(mol.get_bonds_in_molecule(self._bonds_to_break))

        if self._precache_predictions:
            for mol in cache_molecules:
                cache_bonds_to_break.append(
                    mol.get_bonds_in_molecule(self._bonds_to_break)
                )

        self._update_cache(
            molecules + cache_molecules,
            bonds_to_break=bonds_to_break + cache_bonds_to_break,
        )

        for mol, bonds in zip(molecules, bonds_to_break):
            for bond in bonds:
                cache_key = self._cache_key(mol, bond)
                try:
                    predicted_smiles, probabilities, mapped_prod_smiles = self._cache[
                        cache_key
                    ]
                except KeyError:
                    continue

                (
                    predicted_actions,
                    prediction_priors,
                ) = self._actions_from_model_output(
                    mol, predicted_smiles, probabilities, mapped_prod_smiles, bond
                )

                possible_actions.extend(predicted_actions)
                priors.extend(prediction_priors)

        return possible_actions, priors  # type: ignore

    def _actions_from_model_output(
        self,
        mol: TreeMolecule,
        predicted_smiles: List[str],
        probabilities: List[float],
        mapped_product_smiles: Optional[str] = None,
        current_bond: Optional[List[int]] = None,
    ) -> Sequence[SmilesBasedRetroReaction]:
        """
        Create actions given input molecule and output predictions from Chemformer.
        :param mol: input molecules
        :param predicted_smiles: predictions generated by Chemformer
        :param probabilities: prediction probabilities (softmaxed log-likelihoods)
        :param mapped_prod_smiles: the product SMILES mapped with the new atom-mapping
            (used to propagate the old mapping).
        :param current_bond: the bond that was tagged in the input molecule
        :return: list of reactions (actions) and corresponding priors
        """

        actions = []
        priors = []
        for idx, (reactant_smiles, product_smiles) in enumerate(
            zip(predicted_smiles, mapped_product_smiles)
        ):

            metadata = {}
            metadata["policy_probability"] = float(probabilities[idx])
            metadata["policy_probability_rank"] = idx
            metadata["policy_name"] = self.key
            metadata["current_bond"] = current_bond

            reaction = SmilesBasedRetroReaction(
                mol,
                metadata=metadata,
                reactants_str=reactant_smiles,
                mapped_prod_smiles=product_smiles,
            )

            if current_bond not in self._disconnected_bonds(reaction):
                continue

            actions.append(reaction)
            priors.append(float(probabilities[idx]))
        return actions, priors

    def _cache_key(self, molecule: TreeMolecule, bond: Optional[List[int]] = None):
        return f"{molecule.smiles}_{bond}"

    def _make_cache_item(self, prediction: Dict[str, Any]) -> Tuple[Any, ...]:
        return (
            prediction["output"],
            softmax(prediction["lhs"]),
            prediction.get("product_new_mapping"),
        )

    def _make_model_input(
        self,
        molecules: Sequence[TreeMolecule],
        bonds_to_break: Sequence[Sequence[Sequence[int]]],
    ) -> Tuple[List[str], Union[Dict[str, Any], List[str]]]:
        """
        Construct input for the disconnection-aware Chemformer.
        Skip molecules which do not contain bonds to disconnect.
        :param molecules: a list of molecules
        :param bonds_to_break: for each molecule, a list of bonds to break
        :return: model input as a dictionary or None if no molecules to predict on
        """

        if not bonds_to_break:
            return

        product_cache_keys = []
        input_smiles = []
        input_bonds = []
        for molecule, bonds in zip(molecules, bonds_to_break):

            if not bonds:
                continue

            for bond in bonds:
                cache_key = self._cache_key(molecule, bond)

                if cache_key in product_cache_keys or cache_key in self._cache:
                    continue

                input_smiles.append(molecule.mapped_smiles)
                input_bonds.append(bond)
                product_cache_keys.append(cache_key)

        if not product_cache_keys:
            return [], {}

        model_input = {"smiles_list": input_smiles, "bonds_list": input_bonds}
        return product_cache_keys, model_input


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
        input_smiles = []

        for molecule in molecules:
            if molecule.inchi_key in self._cache or molecule.inchi_key in pred_inchis:
                continue
            input_smiles.append(molecule.smiles)
            pred_inchis.append(molecule.inchi_key)

        if not pred_inchis:
            return

        pred_reactants, pred_priors = (
            np.array(item) for item in self.model_zoo.model_call(input_smiles)
        )

        for reactants, priors, inchi in zip(pred_reactants, pred_priors, pred_inchis):
            self._cache[inchi] = (reactants, priors)
