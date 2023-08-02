import glob
import json
import os
import sys
from typing import Dict, List

import pandas as pd
import pytest
import yaml

from aizynthfinder.analysis import RouteCollection
from aizynthfinder.chem import MoleculeException
from aizynthfinder.interfaces import AiZynthApp
from aizynthfinder.interfaces.aizynthapp import main as app_main
from aizynthfinder.interfaces.aizynthcli import main as cli_main
from aizynthfinder.reactiontree import ReactionTree
from aizynthfinder.tools.cat_output import main as cat_main
from aizynthfinder.tools.download_public_data import main as download_main
from aizynthfinder.tools.make_stock import main as make_stock_main

try:
    from aizynthfinder.interfaces.gui import ClusteringGui
except ImportError:
    SUPPORT_CLUSTERING = False
else:
    SUPPORT_CLUSTERING = True


def test_create_gui_app(mocker):
    display_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.display")
    AiZynthApp(configfile="", setup=False)

    display_patch.assert_not_called()

    AiZynthApp(configfile="")

    display_patch.assert_called()


@pytest.mark.xfail(
    condition=not SUPPORT_CLUSTERING, reason="route_distances not installed"
)
def test_create_clustering_gui(mocker, load_reaction_tree):
    collection = RouteCollection(
        reaction_trees=[
            ReactionTree.from_dict(
                load_reaction_tree("routes_for_clustering.json", idx)
            )
            for idx in range(3)
        ]
    )
    display_patch = mocker.patch("aizynthfinder.interfaces.gui.clustering.display")
    ClusteringGui(collection)

    display_patch.assert_called()


@pytest.fixture
def expected_checkpoint_output() -> List[Dict]:
    checkpoint = [
        {
            "processed_smiles": "c1ccccc1",
            "results": {
                "a": 1,
                "b": 2,
                "is_solved": True,
                "stock_info": 1,
                "top_scores": "",
                "trees": 3,
            },
        },
        {
            "processed_smiles": "Cc1ccccc1",
            "results": {
                "a": 1,
                "b": 2,
                "is_solved": True,
                "stock_info": 1,
                "top_scores": "",
                "trees": 3,
            },
        },
        {
            "processed_smiles": "c1ccccc1",
            "results": {
                "a": 1,
                "b": 2,
                "is_solved": True,
                "stock_info": 1,
                "top_scores": "",
                "trees": 3,
            },
        },
        {
            "processed_smiles": "CCO",
            "results": {
                "a": 1,
                "b": 2,
                "is_solved": True,
                "stock_info": 1,
                "top_scores": "",
                "trees": 3,
            },
        },
    ]
    return checkpoint


@pytest.fixture
def multi_smiles_with_checkpoint_results() -> pd.DataFrame:
    results = pd.DataFrame(
        {
            "a": [1, 1, 1, 1],
            "b": [2, 2, 2, 2],
            "is_solved": [True, True, True, True],
            "stock_info": [1, 1, 1, 1],
            "top_scores": ["", "", "", ""],
            "trees": [3, 3, 3, 3],
        },
    )
    return results


def test_app_main_no_output(mocker, tmpdir, add_cli_arguments):
    output_name = str(tmpdir / "temp.ipynb")
    mkstemp_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.tempfile.mkstemp")
    mkstemp_patch.return_value = "dummy", output_name
    popen_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.subprocess.Popen")
    add_cli_arguments("--config config_local.yml")

    app_main()

    popen_patch.assert_called_once()
    with open(output_name, "r") as fileobj:
        lines = fileobj.read()
    assert "from aizynthfinder.interfaces import AiZynthApp" in lines
    assert "config_local.yml" in lines
    assert lines.count("AiZynthApp") == 2


def test_app_main_output(mocker, tmpdir, add_cli_arguments):
    output_name = str(tmpdir / "temp.ipynb")
    popen_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.subprocess.Popen")
    add_cli_arguments("--config config_local.yml --output " + output_name)

    app_main()

    popen_patch.assert_not_called()


def test_cli_single_smiles(mocker, add_cli_arguments, tmpdir, capsys):
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {"a": 1, "b": 2}
    json_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.json.dump")
    output_name = str(tmpdir / "trees.json")
    add_cli_arguments("--smiles COO --config config_local.yml --output " + output_name)

    cli_main()

    finder_patch.assert_called_with(configfile="config_local.yml")
    finder_patch.return_value.expansion_policy.select.assert_called_once()
    finder_patch.return_value.filter_policy.select_all.assert_called_once()
    json_patch.assert_called_once()
    output = capsys.readouterr()
    assert f"Trees saved to {output_name}" in output.out
    assert "Scores for best routes" in output.out
    assert "a: 1" in output.out
    assert "b: 2" in output.out


def test_cli_single_smiles_unsanitizable(mocker, add_cli_arguments, tmpdir, capsys):
    mocker.patch("aizynthfinder.context.config.Configuration.from_file")
    output_name = str(tmpdir / "trees.json")
    add_cli_arguments(
        "--smiles n1(c(=O)[nH]c2c(c1=O)c(c(s2)C(=O)N(C)C)C)c1c(N(=O)O)cccc1 --config config_local.yml --output "
        + output_name
    )
    cli_main()

    output = capsys.readouterr()
    assert "Failed to setup search" in output.out


def test_cli_multiple_smiles(
    mocker,
    add_cli_arguments,
    tmpdir,
    shared_datadir,
    capsys,
    create_dummy_smiles_source,
):
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {
        "a": 1,
        "b": 2,
        "is_solved": True,
    }
    finder_patch.return_value.tree_search.return_value = 1.5
    pd_patch = mocker.patch(
        "aizynthfinder.interfaces.aizynthcli.pd.DataFrame.from_dict"
    )
    smiles_input = create_dummy_smiles_source("txt")
    output_name = str(tmpdir / "data.json.gz")
    add_cli_arguments(
        f"--smiles {smiles_input} --config config_local.yml --output {output_name}"
    )

    cli_main()

    finder_patch.assert_called_with(configfile="config_local.yml")
    pd_patch.assert_called_once()
    output = capsys.readouterr()
    assert output.out.count("Done with") == 4
    assert f"Output saved to {output_name}" in output.out


def test_cli_multiple_smiles_with_empty_checkpoint(
    mocker,
    add_cli_arguments,
    tmpdir,
    shared_datadir,
    create_dummy_smiles_source,
    expected_checkpoint_output,
    multi_smiles_with_checkpoint_results,
):
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {
        "a": 1,
        "b": 2,
        "is_solved": True,
    }
    finder_patch.return_value.tree_search.return_value = 1.5
    finder_patch.return_value.stock_info.return_value = 1
    finder_patch.return_value.trees.return_value = 2
    finder_patch.return_value.routes.dict_with_extra.return_value = 3

    smiles_input = create_dummy_smiles_source("txt")
    output_name = str(tmpdir / "data.json.gz")
    checkpoint = str(tmpdir / "checkpoint.json.gz")
    add_cli_arguments(
        f"--smiles {smiles_input} --config config_local.yml "
        f"--output {output_name} --checkpoint {checkpoint}"
    )

    cli_main()

    with open(checkpoint) as json_file:
        checkpoint_output = [json.loads(line) for line in json_file]
    results_output = pd.read_json(output_name, orient="table")

    assert checkpoint_output == expected_checkpoint_output
    pd.testing.assert_frame_equal(results_output, multi_smiles_with_checkpoint_results)


def test_cli_multiple_smiles_with_checkpoint(
    mocker,
    add_cli_arguments,
    tmpdir,
    shared_datadir,
    create_dummy_smiles_source,
    expected_checkpoint_output,
    multi_smiles_with_checkpoint_results,
):
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {
        "a": 1,
        "b": 2,
        "is_solved": True,
    }
    finder_patch.return_value.tree_search.return_value = 1.5
    finder_patch.return_value.stock_info.return_value = 1
    finder_patch.return_value.trees.return_value = 2
    finder_patch.return_value.routes.dict_with_extra.return_value = 3

    smiles_input = create_dummy_smiles_source("txt")
    output_name = str(tmpdir / "data.json.gz")
    checkpoint = str(tmpdir / "checkpoint.json.gz")

    with open(shared_datadir / "input_checkpoint.json.gz", "r") as from_file, open(
        checkpoint, "w"
    ) as to:
        to.write(from_file.read())

    add_cli_arguments(
        f"--smiles {smiles_input} --config config_local.yml "
        f"--output {output_name} --checkpoint {checkpoint}"
    )

    cli_main()

    with open(checkpoint) as json_file:
        checkpoint_output = [json.loads(line) for line in json_file]
    results_output = pd.read_json(output_name, orient="table")

    assert checkpoint_output == expected_checkpoint_output
    pd.testing.assert_frame_equal(results_output, multi_smiles_with_checkpoint_results)


def test_cli_multiple_smiles_unsanitizable(
    mocker,
    add_cli_arguments,
    tmpdir,
    capsys,
):
    mocker.patch("aizynthfinder.context.config.Configuration.from_file")
    smiles_input = str(tmpdir / "smiles_source.txt")
    with open(smiles_input, "w") as fileobj:
        fileobj.write("n1(c(=O)[nH]c2c(c1=O)c(c(s2)C(=O)N(C)C)C)c1c(N(=O)O)cccc1")
    output_name = str(tmpdir / "data.json.gz")
    add_cli_arguments(
        f"--smiles {smiles_input} --config config_local.yml --output {output_name}"
    )

    cli_main()

    output = capsys.readouterr()
    assert "Failed to setup search" in output.out
    assert output.out.count("Done with") == 0
    assert f"Output saved to {output_name}" in output.out
    assert len(pd.read_json(output_name, orient="table")) == 0


def test_cli_single_smile_with_postprocessing(
    mocker, add_cli_arguments, tmpdir, capsys, shared_datadir
):
    module_path = str(shared_datadir)
    sys.path.append(module_path)
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {"a": 1, "b": 2}
    mocker.patch("aizynthfinder.interfaces.aizynthcli.json.dump")
    output_name = str(tmpdir / "trees.json")
    add_cli_arguments(
        "--post_processing post_processing_test --smiles COO --config config_local.yml --output "
        + output_name
    )

    cli_main()

    output = capsys.readouterr()
    assert "quantity: 5" in output.out
    assert "another quantity: 10" in output.out

    sys.path.remove(module_path)


def test_cli_smiles_argument_incorrect(
    add_cli_arguments,
    tmpdir,
    capsys,
):
    smiles_input = str(tmpdir / "smiles_source.txt")
    add_cli_arguments(f"--smiles {smiles_input} --config config_local.yml")

    cli_main()

    output = capsys.readouterr()
    assert "Cannot start retrosynthesis planning." in output.out
    assert smiles_input in output.out


def test_make_stock_from_plain_file(
    create_dummy_smiles_source, tmpdir, add_cli_arguments, default_config
):
    output_name = str(tmpdir / "temp.hdf5")
    filename = create_dummy_smiles_source("txt")
    add_cli_arguments(f"--files {filename} --output {output_name}")

    make_stock_main()

    default_config.stock.load(filename, "stock1")
    default_config.stock.select(["stock1"])
    assert len(default_config.stock) == 3


def test_cat_main(tmpdir, add_cli_arguments, create_dummy_stock1, create_dummy_stock2):
    filename = str(tmpdir / "output.hdf")
    inputs = [create_dummy_stock1("hdf5"), create_dummy_stock2]
    add_cli_arguments(f"--files {inputs[0]} {inputs[1]} --output {filename}")

    cat_main()

    data = pd.read_hdf(filename, "table")
    assert len(data) == 4


def test_cat_main_json(
    tmpdir, add_cli_arguments, create_dummy_stock1, create_dummy_stock2
):
    filename = str(tmpdir / "output.json.gz")
    inputs = [create_dummy_stock1("hdf5"), create_dummy_stock2]
    add_cli_arguments(f"--files {inputs[0]} {inputs[1]} --output {filename}")

    cat_main()

    data = pd.read_json(filename, orient="table")
    assert len(data) == 4


def test_download_public_data(tmpdir, mocker, add_cli_arguments):
    request_mock = mocker.patch("aizynthfinder.tools.download_public_data.requests.get")
    response_mock = request_mock.return_value
    filecontent = [b"abc", b"def"]
    response_mock.__enter__.return_value.iter_content.return_value = filecontent
    add_cli_arguments(str(tmpdir))

    download_main()

    filenames = glob.glob(str(tmpdir / "*.hdf5"))
    assert len(filenames) == 4
    for filename in filenames:
        with open(filename, "r") as fileobj:
            assert fileobj.read() == "abcdef"

    assert os.path.exists(tmpdir / "config.yml")
    with open(tmpdir / "config.yml", "r") as fileobj:
        config = yaml.load(fileobj.read(), Loader=yaml.SafeLoader)
    policies = config.get("policy", {}).get("files", {})
    assert "uspto" in policies
    assert len(policies["uspto"]) == 2
    stocks = config.get("stock", {}).get("files", {})
    assert "zinc" in stocks
