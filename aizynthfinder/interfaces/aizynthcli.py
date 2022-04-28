""" Module containing classes and routines for the CLI
"""
from __future__ import annotations
import argparse
import json
import os
import warnings
import logging
import importlib
import tempfile
import time
from collections import defaultdict
from typing import TYPE_CHECKING

import pandas as pd

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.utils.files import cat_hdf_files, split_file, start_processes
from aizynthfinder.utils.logging import logger, setup_logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import StrDict, Callable, List, Optional

    _PostProcessingJob = Callable[[AiZynthFinder], StrDict]


def _do_clustering(
    finder: AiZynthFinder,
    results: StrDict,
    detailed_results: bool,
    model_path: str = None,
) -> None:
    time0 = time.perf_counter_ns()
    if model_path:
        kwargs = {"distances_model": "lstm", "model_path": model_path}
    else:
        kwargs = {"distances_model": "ted"}
    results["cluster_labels"] = finder.routes.cluster(n_clusters=0, **kwargs)  # type: ignore
    if not detailed_results:
        return

    results["cluster_time"] = (time.perf_counter_ns() - time0) * 1e-9
    results["distance_matrix"] = finder.routes.distance_matrix().tolist()


def _do_post_processing(
    finder: AiZynthFinder, results: StrDict, jobs: List[_PostProcessingJob]
) -> None:
    for job in jobs:
        results.update(job(finder))


def _get_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser("aizynthcli")
    parser.add_argument(
        "--smiles",
        required=True,
        help="the target molecule smiles or the path of a file containing the smiles",
    )
    parser.add_argument(
        "--config", required=True, help="the filename of a configuration file"
    )
    parser.add_argument(
        "--policy",
        nargs="+",
        default=[],
        help="the name of the expansion policy to use",
    )
    parser.add_argument(
        "--filter", nargs="+", default=[], help="the name of the filter to use"
    )
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
    parser.add_argument(
        "--nproc",
        type=int,
        help="if given, the input is split over a number of processes",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        default=False,
        help="if provided, perform automatic clustering",
    )
    parser.add_argument(
        "--route_distance_model",
        help="if provided, calculate route distances for clustering with this ML model",
    )
    parser.add_argument(
        "--post_processing",
        nargs="+",
        help="a number of modules that performs post-processing tasks",
    )
    return parser.parse_args()


def _load_postprocessing_jobs(modules: Optional[List[str]]) -> List[_PostProcessingJob]:
    jobs: List[_PostProcessingJob] = []
    for module_name in modules or []:
        try:
            module = importlib.import_module(module_name)
        except ModuleNotFoundError:
            pass
        else:
            if hasattr(module, "post_processing"):
                print(f"Adding post-processing job from {module_name}")
                jobs.append(getattr(module, "post_processing"))
    return jobs


def _select_stocks(finder: AiZynthFinder, args: argparse.Namespace) -> None:
    stocks = list(args.stocks)
    try:
        module = importlib.import_module("custom_stock")
    except ModuleNotFoundError:
        pass
    else:
        if hasattr(module, "stock"):
            finder.stock.load(module.stock, "custom_stock")  # type: ignore
            stocks.append("custom_stock")
    finder.stock.select(stocks or finder.stock.items)


def _process_single_smiles(
    smiles: str,
    finder: AiZynthFinder,
    output_name: str,
    do_clustering: bool,
    route_distance_model: Optional[str],
    post_processing: List[_PostProcessingJob],
) -> None:
    output_name = output_name or "trees.json"
    finder.target_smiles = smiles
    try:
        finder.prepare_tree()
    except ValueError as err:
        print(f"Failed to setup search due to: '{str(err).lower()}'")
        return
    finder.tree_search(show_progress=True)
    finder.build_routes()

    with open(output_name, "w") as fileobj:
        json.dump(finder.routes.dicts, fileobj, indent=2)
    logger().info(f"Trees saved to {output_name}")

    scores = ", ".join("%.4f" % score for score in finder.routes.scores)
    logger().info(f"Scores for best routes: {scores}")

    stats = finder.extract_statistics()
    if do_clustering:
        _do_clustering(
            finder, stats, detailed_results=False, model_path=route_distance_model
        )
    _do_post_processing(finder, stats, post_processing)
    stats_str = "\n".join(
        f"{key.replace('_', ' ')}: {value}" for key, value in stats.items()
    )
    logger().info(stats_str)


def _process_multi_smiles(
    filename: str,
    finder: AiZynthFinder,
    output_name: str,
    do_clustering: bool,
    route_distance_model: Optional[str],
    post_processing: List[_PostProcessingJob],
) -> None:
    output_name = output_name or "output.hdf5"
    with open(filename, "r") as fileobj:
        smiles = [line.strip() for line in fileobj.readlines()]

    results = defaultdict(list)
    for smi in smiles:
        finder.target_smiles = smi
        try:
            finder.prepare_tree()
        except ValueError as err:
            print(f"Failed to setup search for {smi} due to: '{str(err).lower()}'")
            continue
        search_time = finder.tree_search()
        finder.build_routes()
        stats = finder.extract_statistics()

        solved_str = "is solved" if stats["is_solved"] else "is not solved"
        logger().info(f"Done with {smi} in {search_time:.3} s and {solved_str}")
        if do_clustering:
            _do_clustering(
                finder, stats, detailed_results=True, model_path=route_distance_model
            )
        _do_post_processing(finder, stats, post_processing)
        for key, value in stats.items():
            results[key].append(value)
        results["top_scores"].append(
            ", ".join("%.4f" % score for score in finder.routes.scores)
        )
        results["trees"].append(finder.routes.dicts)

    data = pd.DataFrame.from_dict(results)
    with warnings.catch_warnings():  # This wil suppress a PerformanceWarning
        warnings.simplefilter("ignore")
        data.to_hdf(output_name, key="table", mode="w")
    logger().info(f"Output saved to {output_name}")


def _multiprocess_smiles(args: argparse.Namespace) -> None:
    def create_cmd(index, filename):
        cmd_args = [
            "aizynthcli",
            "--smiles",
            filename,
            "--config",
            args.config,
            "--output",
            hdf_files[index - 1],
        ]
        if args.policy:
            cmd_args.extend(["--policy"] + args.policy)
        if args.filter:
            cmd_args.extend(["--filter"] + args.filter)
        if args.stocks:
            cmd_args.append("--stocks")
            cmd_args.extend(args.stocks)
        if args.cluster:
            cmd_args.append("--cluster")
        if args.route_distance_model:
            cmd_args.extend(["--route_distance_model", args.route_distance_model])
        if args.post_processing:
            cmd_args.extend(["--post_processing"] + args.post_processing)
        return cmd_args

    if not os.path.exists(args.smiles):
        raise ValueError(
            "For multiprocessing execution the --smiles argument needs to be a filename"
        )

    setup_logger(logging.INFO)
    filenames = split_file(args.smiles, args.nproc)
    hdf_files = [tempfile.mktemp(suffix=".hdf") for _ in range(args.nproc)]
    start_processes(filenames, "aizynthcli", create_cmd)

    if not all(os.path.exists(filename) for filename in hdf_files):
        raise FileNotFoundError(
            "Not all output files produced. Please check the individual log files: 'aizynthcli*.log'"
        )
    cat_hdf_files(hdf_files, args.output or "output.hdf5")


def main() -> None:
    """Entry point for the aizynthcli command"""
    args = _get_arguments()
    if args.nproc:
        _multiprocess_smiles(args)
        return

    multi_smiles = os.path.exists(args.smiles)

    file_level_logging = logging.DEBUG if args.log_to_file else None
    setup_logger(logging.INFO, file_level_logging)

    finder = AiZynthFinder(configfile=args.config)
    _select_stocks(finder, args)
    post_processing = _load_postprocessing_jobs(args.post_processing)
    finder.expansion_policy.select(args.policy or finder.expansion_policy.items[0])
    finder.filter_policy.select(args.filter)

    func = _process_multi_smiles if multi_smiles else _process_single_smiles
    func(
        args.smiles,
        finder,
        args.output,
        args.cluster,
        args.route_distance_model,
        post_processing,
    )


if __name__ == "__main__":
    main()
