""" Module containing classes that interfaces neural network policies
"""
import numpy as np
import pandas as pd

from aizynthfinder.chem import RetroReaction
from aizynthfinder.utils.models import load_model
from aizynthfinder.context.collection import ContextCollection


class PolicyException(Exception):
    """ An exception raised by the Policy classes
    """


def _make_fingerprint(obj, model):
    fingerprint = obj.fingerprint(radius=2, nbits=len(model))
    return fingerprint.reshape([1, len(model)])


class ExpansionPolicy(ContextCollection):
    """
    An abstraction of an expansion policy.

    This policy provides actions (templates) that can be applied to a molecule

    :param config: the configuration of the tree search
    :type config: Configuration
    """

    _collection_name = "expansion policy"

    def __init__(self, config):
        super().__init__()
        self._config = config
        self._stock = config.stock

    def __call__(self, molecules):
        return self.get_actions(molecules)

    def get_actions(self, molecules):
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :type molecules: list of Molecule
        :return: the actions and the priors of those actions
        :rtype: tuple (list of RetroReaction, numpy.ndarray)
        """
        possible_actions = []
        priors = []

        for mol in molecules:
            if mol in self._stock:
                continue

            for policy_key in self.selection:
                model = self[policy_key]["model"]
                templates = self[policy_key]["templates"]

                all_transforms_prop = self._predict(mol, model)
                probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)
                possible_moves = templates.iloc[probable_transforms_idx]
                probs = all_transforms_prop[probable_transforms_idx]

                priors.extend(probs)
                for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
                    metadata = dict(move)
                    del metadata[self._config.template_column]
                    metadata["policy_probability"] = float(probs[idx])
                    metadata["policy_name"] = policy_key
                    metadata["template_code"] = move_index
                    possible_actions.append(
                        RetroReaction(
                            mol, move[self._config.template_column], metadata=metadata
                        )
                    )
        return possible_actions, priors

    def load(self, source, templatefile, key):
        """
        Load a policy and associated templates under the given key

        If `source` is a string, it is taken as a path to a filename and the
        policy is loaded as an `LocalKerasModel` object.

        If `source` is not a string, it is taken as a custom object that
        implements the `__len__` and `predict` methods.

        :param source: the source of the policy model
        :type source: str or object
        :param templatefile: the path to a HDF5 file with the templates
        :type templatefile: str
        :param key: the key or label
        :type key: str
        :raises PolicyException: if the length of the model output vector is not same as the number of templates
        """
        self._logger.info(f"Loading expansion policy model from {source} to {key}")
        model = load_model(source, key)

        self._logger.info(f"Loading templates from {templatefile} to {key}")
        templates = pd.read_hdf(templatefile, "table")

        if hasattr(model, "output_size") and len(templates) != model.output_size:
            raise PolicyException(
                f"The number of templates ({len(templates)}) does not agree with the "
                f"output dimensions of the model ({model.output_size})"
            )

        self._items[key] = {"model": model, "templates": templates}

    def load_from_config(self, **config):
        """
        Load one or more expansion policy from a configuration

        The format should be
        key:
            - path_to_model
            - path_to_templates

        :param config: the configuration
        :type config: key value pairs
        """
        for key, policy_spec in config.items():
            modelfile, templatefile = policy_spec
            self.load(modelfile, templatefile, key)

    def _cutoff_predictions(self, predictions):
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        sortidx = np.argsort(predictions)[::-1]
        cumsum = np.cumsum(predictions[sortidx])
        if any(cumsum >= self._config.cutoff_cumulative):
            maxidx = np.argmin(cumsum < self._config.cutoff_cumulative)
        else:
            maxidx = len(cumsum)
        maxidx = min(maxidx, self._config.cutoff_number) or 1
        return sortidx[:maxidx]

    def _predict(self, mol, model):
        fp_arr = _make_fingerprint(mol, model)
        return np.array(model.predict(fp_arr)).flatten()


class FilterPolicy(ContextCollection):
    """
    An abstraction of a filter policy.

    This policy provides a query on a reaction to determine whether it is feasible

    :param config: the configuration of the tree search
    :type config: Configuration
    """

    _single_selection = True
    _collection_name = "filter policy"

    def __init__(self, config):
        super().__init__()
        self._config = config

    def __call__(self, reaction, return_prob=False):
        return self.is_feasible(reaction, return_prob)

    def is_feasible(self, reaction, return_prob=False):
        """
        Computes if a given reaction is feasible by given
        the reaction fingerprint to a network model

        :param reaction: the reaction to query
        :type reaction: RetroReaction
        :param return_prob: if True, returns both the feasibility and the probability
        :type return_prob: bool, optional
        :return: if the reaction is feasible
        :rtype: bool
        """
        if not self._selection:
            raise PolicyException("No filter policy selected!")

        if not reaction.reactants:
            return False
        prob = self._predict(reaction)
        feasible = prob >= self._config.filter_cutoff
        if return_prob:
            return feasible, prob
        else:
            return feasible

    def load(self, source, key):
        """
        Load a policy under the given key

        If `source` is a string, it is taken as a path to a filename and the
        policy is loaded as an `LocalKerasModel` object.

        If `source` is not a string, it is taken as a custom object that
        implements the `__len__` and `predict` methods.

        :param source: the source of the policy model
        :type source: str or object

        :param key: the key or label
        :type key: str
        """
        self._logger.info(f"Loading filter policy model from {source} to {key}")
        self._items[key] = {"model": load_model(source, key)}

    def load_from_config(self, **config):
        """
        Load one or more filter policy from a configuration

        The format should be:
            key: path_to_model

        :param config: the configuration
        :type config: key value pairs
        """
        for key, filename in config.items():
            self.load(filename, key)

    def _reaction_to_fingerprint(self, reaction, model):
        rxn_fp = _make_fingerprint(reaction, model)
        prod_fp = _make_fingerprint(reaction.mol, model)
        return prod_fp, rxn_fp

    def _predict(self, mol):
        model = self[self.selection]["model"]
        prod_fp, rxn_fp = self._reaction_to_fingerprint(mol, model)
        return model.predict([prod_fp, rxn_fp])[0][0]
