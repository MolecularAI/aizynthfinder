""" Module containing classes that interfaces neural network policies
"""
import numpy as np
import pandas as pd
from deprecated import deprecated

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
        self._selected_policies = []
        self._stock = config.stock

    def __getitem__(self, key):
        return self._policies[key]

    @property
    def selected_policies(self):
        """
        Returns the keys of the selected policy models.

        :return: the selected policy models
        :rtype: str
        """
        return self._selected_policies

    @property
    @deprecated(
        reason="This is being superseded by the selected_policies property",
        version="1.1.0",
    )
    def selected_policy(self):
        """
        Returns the key of the first selected policy model.
        If set it will update the selection of policy.

        :return: the selected policy model
        :rtype: str
        """
        return self._selected_policies[0]

    @selected_policy.setter
    @deprecated(
        reason="This is being removed in favour of the select_policies method",
        version="1.1.0",
    )
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
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

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

            for policy_key in self._selected_policies:
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
                        Reaction(
                            mol, move[self._config.template_column], metadata=metadata
                        )
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
        :raises PolicyException: if the length of the model output vector is not same as the number of templates
        """
        self._logger.info(
            f"Loading policy model from {source} and templates from {templatefile} to {key}"
        )
        model = self._load_model(source)
        templates = pd.read_hdf(templatefile, "table")

        if hasattr(model, "output_size") and len(templates) != model.output_size:
            raise PolicyException(
                f"The number of templates ({len(templates)}) does not agree with the "
                f"output dimensions of the model ({model.output_size})"
            )

        self._policies[key] = {"model": model, "templates": templates}

    def select_policies(self, keys):
        """
        Select the policies to use.

        :param key: the key of the policies to select
        :type key: str or list of str
        """
        if isinstance(keys, str):
            keys = [keys]

        for key in keys:
            if key not in self._policies:
                raise PolicyException(
                    f"Invalid key specified {key} when selecting policy"
                )

        self._selected_policies = list(keys)
        self._logger.info(f"Selected policies: {', '.join(self._selected_policies)}")

    def select_policy(self, key):
        """
        Select another policy to use.

        Use ``select_policies`` to set multiple policies at once.

        :param key: the key of the policy
        :type key: str
        """
        if key not in self._policies:
            raise PolicyException(f"Invalid key specified {key} when selecting policy")

        self._selected_policies.append(key)
        self._logger.info(f"Selected policies: {', '.join(self._selected_policies)}")

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

    def _predict(self, mol, model):
        fp_arr = self._mol_to_fingerprint(mol, model)
        return np.array(model.predict(fp_arr)).flatten()

    @staticmethod
    def _mol_to_fingerprint(mol, model):
        fingerprint = mol.fingerprint(radius=2, nbits=len(model))
        return fingerprint.reshape([1, len(model)])
