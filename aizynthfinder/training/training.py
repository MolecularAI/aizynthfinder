""" Module containing routines to setup the training of policies.
"""
import argparse

from aizynthfinder.training.utils import Config
from aizynthfinder.training.keras_models import train_rollout_keras_model


def main():
    """ Entry-point for the aizynth_training tool
    """
    parser = argparse.ArgumentParser("Tool to train a network policy")
    parser.add_argument("config", help="the filename to a configuration file")
    args = parser.parse_args()

    config = Config(args.config)
    train_rollout_keras_model(config)


if __name__ == "__main__":
    main()
