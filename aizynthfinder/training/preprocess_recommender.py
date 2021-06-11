""" Module routines for pre-processing data for recommender training
"""
import argparse

import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelBinarizer
from scipy import sparse

from aizynthfinder.training.utils import (
    Config,
    split_and_save_data,
    reactants_to_fingerprint,
)


def _get_config() -> Config:
    parser = argparse.ArgumentParser(
        "Tool to pre-process a template library to be used to train a recommender network"
    )
    parser.add_argument("config", help="the filename to a configuration file")
    args = parser.parse_args()

    return Config(args.config)


def _save_unique_templates(dataset: pd.DataFrame, config: Config) -> None:
    dataset = dataset[["retro_template", "template_code"]]
    dataset = dataset.drop_duplicates(subset="template_code", keep="first")
    dataset.set_index("template_code", inplace=True)
    dataset = dataset.sort_index()
    dataset.to_hdf(config.filename("unique_templates"), "table")


def main() -> None:
    """Entry-point for the preprocess_recommender tool"""
    config = _get_config()

    filename = config.filename("library")
    dataset = pd.read_csv(
        filename,
        index_col=False,
        header=None,
        names=config["library_headers"],
    )

    print("Dataset loaded, generating Labels...", flush=True)
    labelbin = LabelBinarizer(neg_label=0, pos_label=1, sparse_output=True)
    labels = labelbin.fit_transform(dataset["template_hash"])
    split_and_save_data(labels, "labels", config)

    print("Labels created and split, generating Inputs...", flush=True)
    reactants = dataset["reactants"].to_numpy()
    inputs = np.apply_along_axis(reactants_to_fingerprint, 0, [reactants], config)
    inputs = sparse.lil_matrix(inputs.T).tocsr()
    split_and_save_data(inputs, "inputs", config)

    print("Inputs created and split, splitting Full Dataset...", flush=True)
    split_and_save_data(dataset, "library", config)

    print("Full Dataset split, creating unique template set", flush=True)
    _save_unique_templates(dataset, config)


if __name__ == "__main__":
    main()
