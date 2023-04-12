""" Module routines for pre-processing data for recommender training
"""
import argparse
from typing import Sequence, Optional

import pandas as pd
import numpy as np

try:
    from sklearn.preprocessing import LabelBinarizer
    from scipy import sparse
except ImportError:
    raise ImportError(
        "Training is not supported by this installation."
        " Please install aizynthfinder with extras dependencies."
    )

from aizynthfinder.training.utils import (
    Config,
    split_and_save_data,
    reactants_to_fingerprint,
    split_reaction_smiles,
)


def _get_config(optional_args: Optional[Sequence[str]] = None) -> Config:
    parser = argparse.ArgumentParser(
        "Tool to pre-process a template library to be used to train a recommender network"
    )
    parser.add_argument("config", help="the filename to a configuration file")
    args = parser.parse_args(optional_args)

    return Config(args.config)


def _save_unique_templates(dataset: pd.DataFrame, config: Config) -> None:
    dataset = dataset[[config["column_map"]["retro_template"], "template_code"]]
    dataset = dataset.drop_duplicates(subset="template_code", keep="first")
    dataset.set_index("template_code", inplace=True)
    dataset = dataset.sort_index()
    dataset.to_hdf(config.filename("unique_templates"), "table")


def main(optional_args: Optional[Sequence[str]] = None) -> None:
    """Entry-point for the preprocess_recommender tool"""
    config = _get_config(optional_args)

    filename = config.filename("library")
    dataset = pd.read_csv(
        filename,
        index_col=False,
        header=0 if config["in_csv_headers"] else None,
        names=None if config["in_csv_headers"] else config["library_headers"],
        sep=config["csv_sep"],
    )
    if config["reaction_smiles_column"]:
        dataset = split_reaction_smiles(dataset, config)

    print("Dataset loaded, generating Labels...", flush=True)
    labelbin = LabelBinarizer(neg_label=0, pos_label=1, sparse_output=True)
    labels = labelbin.fit_transform(dataset[config["column_map"]["template_hash"]])
    split_and_save_data(labels, "labels", config)

    print("Labels created and split, generating Inputs...", flush=True)
    reactants = dataset[config["column_map"]["reactants"]].to_numpy()
    inputs = np.apply_along_axis(reactants_to_fingerprint, 0, [reactants], config)
    inputs = sparse.lil_matrix(inputs.T).tocsr()
    split_and_save_data(inputs, "inputs", config)

    print("Inputs created and split, splitting Full Dataset...", flush=True)
    split_and_save_data(dataset, "library", config)

    print("Full Dataset split, creating unique template set", flush=True)
    _save_unique_templates(dataset, config)


if __name__ == "__main__":
    main()
