""" Module containing classes and routines for the CLI
"""
import argparse
import json
import os
import warnings
import logging
import importlib
from collections import defaultdict

import pandas as pd

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.utils.logging import logger, setup_logger


def _get_arguments():
    parser = argparse.ArgumentParser("aizynthcli")
    parser.add_argument(
        "--smiles",
        required=True,
        help="the target molecule smiles or the path of a file containing the smiles",
    )
    parser.add_argument(
        "--config", required=True, help="the filename of a configuration file"
    )
    parser.add_argument("--policy", default="", help="the name of the policy to use")
    parser.add_argument(
        "--stocks", nargs="+", default=[], help="the name of the stocks to use"
    )
    parser.add_argument(
        "--output", help="the name of the output file (JSON or HDF5 file)"
    )
    parser.add_argument(
        "--log_to_file",
        action="store_true",
        default=False,
        help="if provided, detailed logging to file is enabled",
    )
    return parser.parse_args()


def _select_stocks(finder, args):
    stocks = list(args.stocks)
    try:
        module = importlib.import_module("custom_stock")
    except ModuleNotFoundError:
        pass
    else:
        if hasattr(module, "stock"):
            finder.stock.load_stock(module.stock, "custom_stock")
            stocks.append("custom_stock")
    finder.stock.select_stocks(stocks or finder.stock.available_stocks())


def _process_single_smiles(smiles, finder, output_name):
    output_name = output_name or "trees.json"
    finder.target_smiles = smiles
    finder.prepare_tree()
    finder.tree_search(show_progress=True)
    finder.build_routes()

    with open(output_name, "w") as fileobj:
        json.dump(finder.routes.dicts, fileobj, indent=2)
    logger().info(f"Trees saved to {output_name}")

    scores = ", ".join("%.4f" % score for score in finder.routes.scores)
    logger().info(f"Scores for best routes: {scores}")

    stats = finder.extract_statistics()
    stats_str = "\n".join(
        f"{key.replace('_', ' ')}: {value}" for key, value in stats.items()
    )
    logger().info(stats_str)


def _process_multi_smiles(filename, finder, output_name):
    output_name = output_name or "output.hdf5"
    with open(filename, "r") as fileobj:
        smiles = [line.strip() for line in fileobj.readlines()]

    results = defaultdict(list)
    for smi in smiles:
        finder.target_smiles = smi
        finder.prepare_tree()
        search_time = finder.tree_search()
        finder.build_routes()
        stats = finder.extract_statistics()

        logger().info(f"Done with {smi} in {search_time:.3} s")
        for key, value in stats.items():
            results[key].append(value)
        results["top_scores"].append(
            ", ".join("%.4f" % score for score in finder.routes.scores)
        )
        results["trees"].append(finder.routes.dicts)

    data = pd.DataFrame.from_dict(results)
    with warnings.catch_warnings():  # This wil supress a PerformanceWarning
        warnings.simplefilter("ignore")
        data.to_hdf(output_name, key="table", mode="w")
    logger().info(f"Output saved to {output_name}")


def main():
    """ Entry point for the aizynthcli command
    """
    args = _get_arguments()
    multi_smiles = os.path.exists(args.smiles)

    file_level_logging = logging.DEBUG if args.log_to_file else None
    setup_logger(logging.INFO, file_level_logging)

    finder = AiZynthFinder(configfile=args.config)
    _select_stocks(finder, args)
    finder.policy.select_policy(args.policy or finder.policy.available_policies()[0])

    if multi_smiles:
        _process_multi_smiles(args.smiles, finder, args.output)
    else:
        _process_single_smiles(args.smiles, finder, args.output)


if __name__ == "__main__":
    main()
