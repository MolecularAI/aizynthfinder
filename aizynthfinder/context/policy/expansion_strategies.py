""" Module containing classes that implements different expansion policy strategies
"""
from __future__ import annotations
import abc
from re import template
from typing import TYPE_CHECKING

import os
import numpy as np
import pandas as pd
import gcsfs
import h5py

from aizynthfinder.chem import TemplatedRetroReaction
from aizynthfinder.utils.models import load_model
from aizynthfinder.utils.logging import logger
from aizynthfinder.context.policy.utils import _make_fingerprint
from aizynthfinder.utils.exceptions import PolicyException

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Any, Sequence, List, Tuple
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction

PROJECT_NAME = os.environ.get("PROJECT_NAME")
GOOGLE_APPLICATION_CREDENTIALS = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS")

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
        self.inchi_fail = 0
        self.rdkit_fail = 0

    def __call__(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[RetroReaction], List[float]]:
        return self.get_actions(molecules)

    @abc.abstractmethod
    def get_actions(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules

        :param molecules: the molecules to consider
        :return: the actions and the priors of those actions
        """


class TemplateBasedExpansionStrategy(ExpansionStrategy):
    """
    A template-based expansion strategy that will return `TemplatedRetroReaction` objects upon expansion.

    :param key: the key or label
    :param config: the configuration of the tree search
    :param source: the source of the policy model
    :param templatefile: the path to a HDF5 file with the templates
    :raises PolicyException: if the length of the model output vector is not same as the number of templates
    """

    _required_kwargs = ["source", "templatefile"]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        super().__init__(key, config, **kwargs)

        source = kwargs["source"]
        templatefile = kwargs["templatefile"]
        # print(f"Source = {source}\tTemplate = {templatefile}")

        self._logger.info(
            f"Loading template-based expansion policy model from {source} to {self.key}"
        )
        # print(f"\n\nFrom TemplateBasedExpansionStrategy:\nsource = {source}\nkey = {self.key}\nconfig = {self._config}")
        self.model = load_model(source, self.key, self._config.use_remote_models)

        self._logger.info(f"Loading templates from {templatefile} to {self.key}")
        print(f"Loading templates from {templatefile} to {self.key} and {source}")
        ext = templatefile.rsplit(".")[-1]
        if "gs://" in templatefile:
            # Load from google cloud
            bucket = gcsfs.GCSFileSystem(
                token=GOOGLE_APPLICATION_CREDENTIALS,
                project=PROJECT_NAME)
            with bucket.open(templatefile, "rb") as cloud_file:
                if ext == "hdf5":
                    # Soln from: https://stackoverflow.com/questions/40472912/hdf5-file-to-pandas-dataframe
                    # but doesn't seem to work properly
                    data = h5py.File(cloud_file, "r")
                    try:
                        self.templates: pd.DataFrame = pd.DataFrame(np.array(data["table"]))
                    except KeyError:
                        pass
                        # raise KeyError(f"Error with {key} and {source}.\nFile keys = {data.keys()}")
                elif ext == "csv":
                    self.templates: pd.DataFrame = pd.read_csv(cloud_file, index_col=0)
        elif ext == "csv":
            self.templates: pd.DataFrame = pd.read_csv(templatefile)
        else:       
            self.templates: pd.DataFrame = pd.read_hdf(templatefile, "table")

        if hasattr(self.model, "output_size") and len(self.templates) != self.model.output_size:  # type: ignore
            raise PolicyException(
                f"The number of templates ({len(self.templates)}) does not agree with the "  # type: ignore
                f"output dimensions of the model ({self.model.output_size})"
            )

    # pylint: disable=R0914
    def get_actions(
        self, molecules: Sequence[TreeMolecule]
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :return: the actions and the priors of those actions
        """

        possible_actions = []
        priors = []

        for mol in molecules:
            model = self.model
            templates = self.templates

            all_transforms_prop = self._predict(mol, model)
            probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)
            possible_moves = templates.iloc[probable_transforms_idx]
            probs = all_transforms_prop[probable_transforms_idx]

            priors.extend(probs)
            for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
                metadata = dict(move)
                del metadata[self._config.template_column]
                metadata["policy_probability"] = float(probs[idx].round(4))
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key
                metadata["template_code"] = move_index
                metadata["template"] = move[self._config.template_column]
                possible_actions.append(
                    TemplatedRetroReaction(
                        mol,
                        smarts=move[self._config.template_column],
                        metadata=metadata,
                        reaction_source=self._config.reaction_source,
                        template_fallback=self._config.template_fallback,
                        use_rdchiral=self._config.use_rdchiral,
                        templates=self.templates,
                        inchi_fail=self.inchi_fail,
                        rdkit_fail=self.rdkit_fail,
                    )
                )
        return possible_actions, priors  # type: ignore

    def _cutoff_predictions(self, predictions: np.ndarray) -> np.ndarray:
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        sortidx = np.argsort(predictions)[::-1]
        cumsum: np.ndarray = np.cumsum(predictions[sortidx])
        if any(cumsum >= self._config.cutoff_cumulative):
            maxidx = int(np.argmin(cumsum < self._config.cutoff_cumulative))
        else:
            maxidx = len(cumsum)
        maxidx = min(maxidx, self._config.cutoff_number) or 1
        return sortidx[:maxidx]

    @staticmethod
    def _predict(mol: TreeMolecule, model: Any) -> np.ndarray:
        fp_arr = _make_fingerprint(mol, model)
        return np.array(model.predict(fp_arr)).flatten()
