import os
import gzip
import json

import pytest
import pandas as pd

from aizynthfinder.utils.files import (
    cat_datafiles,
    split_file,
    start_processes,
    read_datafile,
    save_datafile,
)


@pytest.fixture
def create_dummy_file(tmpdir, mocker):
    patched_tempfile = mocker.patch("aizynthfinder.utils.files.tempfile.mktemp")
    split_files = [
        tmpdir / "split1",
        tmpdir / "split2",
        tmpdir / "split3",
    ]
    patched_tempfile.side_effect = split_files
    filename = tmpdir / "input"

    def wrapper(content):
        with open(filename, "w") as fileobj:
            fileobj.write(content)
        return filename, split_files

    return wrapper


def test_split_file_even(create_dummy_file):
    filename, split_files = create_dummy_file("\n".join(list("abcdef")))

    split_file(filename, 3)

    read_lines = []
    for filename in split_files:
        assert os.path.exists(filename)
        with open(filename, "r") as fileobj:
            read_lines.append(fileobj.read())
    assert read_lines[0] == "a\nb"
    assert read_lines[1] == "c\nd"
    assert read_lines[2] == "e\nf"


def test_split_file_odd(create_dummy_file):
    filename, split_files = create_dummy_file("\n".join(list("abcdefg")))

    split_file(filename, 3)

    read_lines = []
    for filename in split_files:
        assert os.path.exists(filename)
        with open(filename, "r") as fileobj:
            read_lines.append(fileobj.read())
    assert read_lines[0] == "a\nb\nc"
    assert read_lines[1] == "d\ne"
    assert read_lines[2] == "f\ng"


def test_start_processes(tmpdir):

    script_filename = str(tmpdir / "dummy.py")
    with open(script_filename, "w") as fileobj:
        fileobj.write("import sys\nimport time\nprint(sys.argv[1])\ntime.sleep(2)\n")

    def create_cmd(index, filename):
        return ["python", script_filename, f"{filename}-{index}"]

    start_processes(["dummy", "dummy"], str(tmpdir / "log"), create_cmd, 2)

    for index in [1, 2]:
        logfile = str(tmpdir / f"log{index}.log")
        assert os.path.exists(logfile)
        with open(logfile, "r") as fileobj:
            lines = fileobj.read()
        assert lines == f"dummy-{index}\n"


def test_cat_hdf(create_dummy_stock1, create_dummy_stock2, tmpdir):
    filename = str(tmpdir / "output.hdf")
    inputs = [create_dummy_stock1("hdf5"), create_dummy_stock2]

    cat_datafiles(inputs, filename)

    data = pd.read_hdf(filename, "table")
    assert len(data) == 4
    assert list(data.inchi_key.values) == [
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "YXFVVABEGXRONW-UHFFFAOYSA-N",
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "ISWSIDIOOBJBQZ-UHFFFAOYSA-N",
    ]


def test_cat_hdf_no_trees(tmpdir, create_dummy_stock1, create_dummy_stock2):
    hdf_filename = str(tmpdir / "output.hdf")
    tree_filename = str(tmpdir / "trees.json")
    inputs = [create_dummy_stock1("hdf5"), create_dummy_stock2]

    cat_datafiles(inputs, hdf_filename, tree_filename)

    assert not os.path.exists(tree_filename)


def test_cat_hdf_trees(tmpdir):
    hdf_filename = str(tmpdir / "output.hdf")
    tree_filename = str(tmpdir / "trees.json")
    filename1 = str(tmpdir / "file1.hdf5")
    filename2 = str(tmpdir / "file2.hdf5")
    trees1 = [[1], [2]]
    trees2 = [[3], [4]]
    pd.DataFrame({"mol": ["A", "B"], "trees": trees1}).to_hdf(filename1, "table")
    pd.DataFrame({"mol": ["A", "B"], "trees": trees2}).to_hdf(filename2, "table")

    cat_datafiles([filename1, filename2], hdf_filename, tree_filename)

    assert os.path.exists(tree_filename + ".gz")
    with gzip.open(tree_filename + ".gz", "rt", encoding="UTF-8") as fileobj:
        trees_cat = json.load(fileobj)
    assert trees_cat == trees1 + trees2
    assert "trees" not in pd.read_hdf(hdf_filename, "table")


@pytest.mark.parametrize(
    ("filename"),
    [
        ("temp.json"),
        ("temp.hdf5"),
    ],
)
def test_save_load_datafile_roundtrip(filename, tmpdir):
    data1 = pd.DataFrame({"a": [0, 1, 2], "b": [2, 3, 4]})

    save_datafile(data1, tmpdir / filename)

    data2 = read_datafile(tmpdir / filename)

    assert data1.columns.to_list() == data2.columns.to_list()
    assert data1.a.to_list() == data2.a.to_list()
    assert data1.b.to_list() == data2.b.to_list()
