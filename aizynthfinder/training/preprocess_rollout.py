""" Module routines for pre-processing data for rollout policy training
"""
import argparse
import os

import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelBinarizer, LabelEncoder
from scipy import sparse

from aizynthfinder.training.utils import (
    Config,
    split_and_save_data,
    smiles_to_fingerprint,
)


def _filter_dataset(config):

    filename = config.filename("raw_library")
    if not os.path.exists(filename):
        raise FileNotFoundError(
            f"The file {filename} is missing - cannot proceed without the full template library."
        )

    # Skipping the last header as it is not available in the raw data
    full_data = pd.read_csv(
        filename, index_col=False, header=None, names=config["library_headers"][:-1],
    )
    full_data = full_data.drop_duplicates(subset="reaction_hash")
    template_group = full_data.groupby("template_hash")
    template_group = template_group.size().sort_values(ascending=False)
    min_index = template_group[template_group >= config["template_occurance"]].index
    dataset = full_data[full_data["template_hash"].isin(min_index)]

    template_labels = LabelEncoder()
    dataset["template_code"] = template_labels.fit_transform(dataset["template_hash"])
    dataset.to_csv(
        config.filename("library"), mode="w", header=False, index=False,
    )
    return dataset


def _get_config():
    parser = argparse.ArgumentParser(
        "Tool to pre-process a template library to be used in training a rollout network policy"
    )
    parser.add_argument("config", help="the filename to a configuration file")
    args = parser.parse_args()

    return Config(args.config)


def _save_unique_templates(dataset, config):
    template_group = dataset.groupby("template_hash", sort=False).size()
    dataset = dataset[["retro_template", "template_code"] + config["metadata_headers"]]
    dataset = dataset.drop_duplicates(subset="template_code", keep="first")
    dataset["library_occurence"] = template_group.values
    dataset.set_index("template_code", inplace=True)
    dataset = dataset.sort_index()
    dataset.to_hdf(config.filename("unique_templates"), "table")


def main():
    """ Entry-point for the preprocess_rollout tool
    """
    config = _get_config()

    filename = config.filename("library")
    if not os.path.exists(filename):
        dataset = _filter_dataset(config)
    else:
        dataset = pd.read_csv(
            filename, index_col=False, header=None, names=config["library_headers"],
        )

    print("Dataset filtered/loaded, generating labels...", flush=True)
    lb = LabelBinarizer(neg_label=0, pos_label=1, sparse_output=True)
    labels = lb.fit_transform(dataset["template_hash"])
    split_and_save_data(labels, "labels", config)

    print("Labels created and splitted, generating inputs...", flush=True)
    products = dataset["products"].to_numpy()
    inputs = np.apply_along_axis(smiles_to_fingerprint, 0, [products], config)
    inputs = sparse.lil_matrix(inputs.T).tocsr()
    split_and_save_data(inputs, "inputs", config)

    print("Inputs created and splitted, splitting full Dataset...", flush=True)
    split_and_save_data(dataset, "library", config)

    print("Full Dataset splitted, creating unique template set", flush=True)
    _save_unique_templates(dataset, config)


if __name__ == "__main__":
    main()
