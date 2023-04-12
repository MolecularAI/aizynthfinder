""" Module routines for pre-processing data for filter policy training
"""
import argparse
from typing import Sequence, Optional

import pandas as pd
import numpy as np

try:
    from scipy import sparse
except ImportError:
    raise ImportError(
        "Training is not supported by this installation."
        " Please install aizynthfinder with extras dependencies."
    )

from aizynthfinder.training.utils import (
    Config,
    split_and_save_data,
    smiles_to_fingerprint,
    reaction_to_fingerprints,
    split_reaction_smiles,
)


def _get_config(optional_args: Optional[Sequence[str]] = None) -> Config:
    parser = argparse.ArgumentParser(
        "Tool to pre-process a template library to be used to train a in-scope filter network policy"
    )
    parser.add_argument("config", help="the filename to a configuration file")
    args = parser.parse_args(optional_args)

    return Config(args.config)


def main(optional_args: Optional[Sequence[str]] = None) -> None:
    """Entry-point for the preprocess_filter tool"""
    config = _get_config(optional_args)

    true_dataset = pd.read_csv(
        config.filename("library"),
        index_col=False,
        header=0 if config["in_csv_headers"] else None,
        names=None if config["in_csv_headers"] else config["library_headers"][:-1],
        sep=config["csv_sep"],
    )
    true_dataset["true_product"] = 1
    false_dataset = pd.read_csv(
        config.filename("false_library"),
        index_col=False,
        header=0 if config["in_csv_headers"] else None,
        names=None if config["in_csv_headers"] else config["library_headers"][:-1],
        sep=config["csv_sep"],
    )
    false_dataset["true_product"] = 0
    dataset = pd.concat([true_dataset, false_dataset])

    if config["reaction_smiles_column"]:
        dataset = split_reaction_smiles(dataset, config)

    print("Dataset loaded, generating Labels...", flush=True)
    labels = dataset["true_product"].to_numpy()
    split_and_save_data(labels, "labels", config)

    print("Labels created and split, generating Inputs...", flush=True)
    products = dataset[config["column_map"]["products"]].to_numpy()
    reactants = dataset[config["column_map"]["reactants"]].to_numpy()
    inputs = np.apply_along_axis(
        reaction_to_fingerprints, 0, [products, reactants], config
    ).astype(np.int8)
    inputs = sparse.lil_matrix(inputs.T).tocsr()
    split_and_save_data(inputs, "inputs2", config)

    inputs = np.apply_along_axis(smiles_to_fingerprint, 0, [products], config).astype(
        np.int8
    )
    inputs = sparse.lil_matrix(inputs.T).tocsr()
    split_and_save_data(inputs, "inputs", config)

    print("Inputs created and split, splitting Full Dataset...", flush=True)
    split_and_save_data(dataset, "library", config)


if __name__ == "__main__":
    main()
