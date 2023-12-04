"""Module containing routines to work with files and processes."""
from __future__ import annotations

import gzip
import json
import subprocess
import tempfile
import time
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from deprecated import deprecated

from aizynthfinder.utils.logging import logger

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Any,
        Callable,
        List,
        Optional,
        Sequence,
        Union,
    )


def read_datafile(filename: Union[str, Path]) -> pd.DataFrame:
    """
    Read aizynth output from disc in either .hdf5 or .json format

    :param filename: the path to the data
    :return: the loaded data
    """
    filename_str = str(filename)
    if filename_str.endswith(".hdf5") or filename_str.endswith(".hdf"):
        return pd.read_hdf(filename, "table")
    return pd.read_json(filename, orient="table")


def save_datafile(data: pd.DataFrame, filename: Union[str, Path]) -> None:
    """
    Save the given data to disc in either .hdf5 or .json format

    :param data: the data to save
    :param filename: the path to the data
    """
    filename_str = str(filename)
    if filename_str.endswith(".hdf5") or filename_str.endswith(".hdf"):
        with warnings.catch_warnings():  # This wil suppress a PerformanceWarning
            warnings.simplefilter("ignore")
            data.to_hdf(filename, key="table")
    else:
        data.to_json(filename, orient="table")


@deprecated(version="4.0.0", reason="replaced by 'cat_datafiles'")
def cat_hdf_files(
    input_files: List[str], output_name: str, trees_name: Optional[str] = None
) -> None:
    """Concatenate hdf5 or json datafiles"""
    cat_datafiles(input_files, output_name, trees_name)


def cat_datafiles(
    input_files: List[str], output_name: str, trees_name: Optional[str] = None
) -> None:
    """
    Concatenate hdf5 or json datafiles

    if `tree_name` is given, will take out the `trees` column
    from the tables and save it to a gzipped-json file.

    :param input_files: the paths to the files to concatenate
    :param output_name: the name of the concatenated file
    :param trees_name: the name of the concatenated trees
    """
    data = read_datafile(input_files[0])
    if "trees" not in data.columns:
        trees_name = None

    if trees_name:
        columns = list(data.columns)
        columns.remove("trees")
        trees = list(data["trees"].values)
        data = data[columns]

    for filename in input_files[1:]:
        new_data = read_datafile(filename)
        if trees_name:
            trees.extend(new_data["trees"].values)
            new_data = new_data[columns]
        data = pd.concat([data, new_data])

    save_datafile(data.reset_index(drop=True), output_name)
    if trees_name:
        if not trees_name.endswith(".gz"):
            trees_name += ".gz"
        with gzip.open(trees_name, "wt", encoding="UTF-8") as fileobj:
            json.dump(trees, fileobj)


def split_file(filename: str, nparts: int) -> List[str]:
    """
    Split the content of a text file into a given number of temporary files

    :param filename: the path to the file to split
    :param nparts: the number of parts to create
    :return: list of filenames of the parts
    """
    with open(filename, "r") as fileobj:
        lines = fileobj.read().splitlines()

    filenames = []
    batch_size, remainder = divmod(len(lines), nparts)
    stop = 0
    for part in range(1, nparts + 1):
        start = stop
        stop += batch_size + 1 if part <= remainder else batch_size
        filenames.append(tempfile.mktemp())
        with open(filenames[-1], "w") as fileobj:
            fileobj.write("\n".join(lines[start:stop]))
    return filenames


def start_processes(
    inputs: Sequence[Any], log_prefix: str, cmd_callback: Callable, poll_freq: int = 5
) -> None:
    """
    Start a number of background processes and wait for them
    to complete.

    The standard output and standard error is saved to a log file.

    The command to start for each process is given by the ``cmd_callback``
    function that takes the index of the process and an item of the input
    as arguments.

    :param inputs: a sequence of input to the processes
    :param log_prefix: the prefix to the log file of each processes
    :param cmd_callback: function that creates the process commands
    :param poll_freq: the polling frequency for checking if processes are completed
    """
    processes = []
    output_fileobjs = []
    for index, iinput in enumerate(inputs, 1):
        output_fileobjs.append(open(f"{log_prefix}{index}.log", "w"))
        cmd = cmd_callback(index, iinput)
        processes.append(
            subprocess.Popen(cmd, stdout=output_fileobjs[-1], stderr=subprocess.STDOUT)
        )
        logger().info(f"Started background task with pid={processes[-1].pid}")

    logger().info("Waiting for background tasks to complete...")
    not_finished = True
    while not_finished:
        time.sleep(poll_freq)
        not_finished = False
        for process, fileobj in zip(processes, output_fileobjs):
            fileobj.flush()
            if process.poll() is None:
                not_finished = True

    for fileobj in output_fileobjs:
        fileobj.close()
