""" Module containing routines to setup the training of policies.
"""
import argparse
from typing import Optional, Sequence

from aizynthfinder.training.utils import Config
from aizynthfinder.training.keras_models import (
    train_expansion_keras_model,
    train_filter_keras_model,
    train_recommender_keras_model,
)


def main(optional_args: Optional[Sequence[str]] = None) -> None:
    """Entry-point for the aizynth_training tool"""
    parser = argparse.ArgumentParser("Tool to train a network policy")
    parser.add_argument("config", help="the filename to a configuration file")
    parser.add_argument(
        "model",
        choices=["expansion", "filter", "recommender"],
        help="the model to train",
    )
    args = parser.parse_args(optional_args)

    config = Config(args.config)
    if args.model == "expansion":
        train_expansion_keras_model(config)
    elif args.model == "filter":
        train_filter_keras_model(config)
    elif args.model == "recommender":
        train_recommender_keras_model(config)


if __name__ == "__main__":
    main()
