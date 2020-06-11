""" Module containing classes that interfaces neural network policies
"""
import numpy as np
import pandas as pd

from aizynthfinder.chem import Reaction
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.keras_utils import LocalKerasModel


class PolicyException(Exception):
    """ An exception raised by the Policy classes
    """


class Policy:
    """
    An abstraction of the network policy that prioritizes
    the tree exploration.

    One can obtain individual policies with:

    .. code-block::

        a_policy = policy["key"]


    Individual policies are dictionaries with the keys
    "model" and "templates", representing the loaded Keras
    model and the loaded template dataset.

    :param config: settings of the tree search algorithm
    :type config: Configuration
    """

    def __init__(self, config):
        self._config = config
        self._logger = logger()
        self._policies = {}
        self._selected_policy = None
        self._policy_model = None
        self._stock = config.stock
        self._templates = None

    def __getitem__(self, key):
        return self._policies[key]

    @property
    def selected_policy(self):
        """
        Returns the key of the selected policy model.
        If set it will update the selection of policy.

        :return: the selected policy model
        :rtype: str
        """
        return self._selected_policy

    @selected_policy.setter
    def selected_policy(self, key):
        self.select_policy(key)

    def available_policies(self):
        """
        Return the keys for the available policies

        :return: the policy keys
        :rtype: tuple of str
        """
        return tuple(self._policies.keys())

    def get_actions(self, molecules):
        """
        Get all the probable actions of a set of molecules, using the policy and given cutoffs

        :param molecules: the molecules to consider
        :type molecules: list of Molecule
        :return: the actions and the priors of those actions
        :rtype: tuple (list of Reaction, numpy.ndarray)
        """
        possible_actions = []
        priors = []

        for mol in molecules:
            if mol in self._stock:
                continue
            all_transforms_prop = self._predict(mol)
            probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)
            possible_moves = self._templates.iloc[probable_transforms_idx]
            probs = all_transforms_prop[probable_transforms_idx]

            priors.extend(probs)
            for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
                metadata = dict(move)
                del metadata[self._config.template_column]
                metadata["policy_probability"] = float(probs[idx])
                metadata["template_code"] = move_index

                possible_actions.append(
                    Reaction(mol, move[self._config.template_column], metadata=metadata)
                )
        return possible_actions, priors

    def load_policy(self, source, templatefile, key):
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
        """
        self._logger.info(
            f"Loading policy model from {source} and templates from {templatefile} to {key}"
        )
        model = self._load_model(source)
        templates = pd.read_hdf(templatefile, "table")
        self._policies[key] = {"model": model, "templates": templates}

    def select_policy(self, key):
        """
        Select and prepare the policy

        :param key: the key of the policy
        :type key: str
        """
        if key not in self._policies:
            raise PolicyException(f"Invalid key specified {key} when selecting policy")

        self._selected_policy = key
        self._policy_model = self._policies[key]["model"]
        self._templates = self._policies[key]["templates"]

        self._logger.info(f"Selected policy is {key}")

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

    def _load_model(self, source):
        if isinstance(source, str):
            return LocalKerasModel(source)
        return source

    def _mol_to_fingerprint(self, mol):
        fingerprint = mol.fingerprint(radius=2, nbits=len(self._policy_model))
        return fingerprint.reshape([1, len(self._policy_model)])

    def _predict(self, mol):
        fp_arr = self._mol_to_fingerprint(mol)
        return np.array(self._policy_model.predict(fp_arr)).flatten()
