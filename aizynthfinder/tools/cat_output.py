""" Module containing a CLI for concatenating output files (hdf5/json.gz files)
"""
import argparse

from aizynthfinder.utils.files import cat_datafiles


def main() -> None:
    """Entry-point for the cat_aizynth_output CLI"""
    parser = argparse.ArgumentParser("cat_aizynthcli_output")
    parser.add_argument(
        "--files",
        required=True,
        nargs="+",
        help="the output filenames",
    )
    parser.add_argument(
        "--output",
        required=True,
        default="output.json.gz",
        help="the name of the concatenate output file ",
    )
    parser.add_argument(
        "--trees",
        help="if given, save all trees to this file",
    )
    args = parser.parse_args()

    cat_datafiles(args.files, args.output, args.trees)


if __name__ == "__main__":
    main()
