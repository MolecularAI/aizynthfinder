import os
import glob

import pandas as pd
import yaml
import pytest

from aizynthfinder.interfaces import AiZynthApp
from aizynthfinder.interfaces.gui import ClusteringGui
from aizynthfinder.interfaces.aizynthapp import main as app_main
from aizynthfinder.interfaces.aizynthcli import main as cli_main
from aizynthfinder.tools.make_stock import main as make_stock_main
from aizynthfinder.tools.cat_output import main as cat_main
from aizynthfinder.training.preprocess_expansion import main as expansion_main
from aizynthfinder.training.preprocess_recommender import main as recommender_main
from aizynthfinder.training.preprocess_filter import main as filter_main
from aizynthfinder.training.make_false_products import main as make_false_main
from aizynthfinder.tools.download_public_data import main as download_main
from aizynthfinder.training.utils import Config
from aizynthfinder.chem import MoleculeException
from aizynthfinder.analysis import RouteCollection, ReactionTree


def test_create_gui_app(mocker):
    display_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.display")
    AiZynthApp(configfile=None, setup=False)

    display_patch.assert_not_called()

    AiZynthApp(configfile=None)

    display_patch.assert_called()


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
    json_patch.assert_called_once()
    output = capsys.readouterr()
    assert f"Trees saved to {output_name}" in output.out
    assert "Scores for best routes" in output.out
    assert "a: 1" in output.out
    assert "b: 2" in output.out


def test_cli_multiple_smiles(mocker, add_cli_arguments, tmpdir, shared_datadir, capsys):
    finder_patch = mocker.patch("aizynthfinder.interfaces.aizynthcli.AiZynthFinder")
    finder_patch.return_value.extract_statistics.return_value = {"a": 1, "b": 2}
    finder_patch.return_value.tree_search.return_value = 1.5
    pd_patch = mocker.patch(
        "aizynthfinder.interfaces.aizynthcli.pd.DataFrame.from_dict"
    )
    smiles_input = str(shared_datadir / "smiles_plain.txt")
    output_name = str(tmpdir / "data.hdf5")
    add_cli_arguments(
        f"--smiles {smiles_input} --config config_local.yml --output {output_name}"
    )

    cli_main()

    finder_patch.assert_called_with(configfile="config_local.yml")
    pd_patch.assert_called_once()
    output = capsys.readouterr()
    assert output.out.count("Done with") == 4
    assert f"Output saved to {output_name}" in output.out


def test_make_stock_from_plain_file(shared_datadir, tmpdir, add_cli_arguments, stock):
    output_name = str(tmpdir / "temp.hdf5")
    filename = str(shared_datadir / "smiles_plain.txt")
    add_cli_arguments(f"--files {filename} --output {output_name}")

    make_stock_main()

    stock.load(filename, "stock1")
    stock.select(["stock1"])
    assert len(stock) == 3


def test_cat_main(shared_datadir, tmpdir, add_cli_arguments):
    filename = str(tmpdir / "output.hdf")
    inputs = [str(shared_datadir / "stock1.h5"), str(shared_datadir / "stock2.h5")]
    add_cli_arguments(f"--files {inputs[0]} {inputs[1]} --output {filename}")

    cat_main()

    data = pd.read_hdf(filename, "table")
    assert len(data) == 4


def test_preprocess_expansion(write_yaml, shared_datadir, add_cli_arguments):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "dummy"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
        }
    )
    add_cli_arguments(config_path)

    expansion_main()

    with open(shared_datadir / "dummy_template_library.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 10

    with open(shared_datadir / "dummy_training.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 6

    with open(shared_datadir / "dummy_testing.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2

    with open(shared_datadir / "dummy_validation.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2

    data = pd.read_hdf(shared_datadir / "dummy_unique_templates.hdf5", "table")
    config = Config(config_path)
    assert len(data) == 2
    assert "retro_template" in data.columns
    assert "library_occurence" in data.columns
    for column in config["metadata_headers"]:
        assert column in data.columns


def test_preprocess_expansion_bad_product(
    write_yaml, shared_datadir, add_cli_arguments
):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "dummy_sani"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
        }
    )
    add_cli_arguments(config_path)
    with pytest.raises(MoleculeException):
        expansion_main()


def test_preprocess_expansion_skip_bad_product(
    write_yaml, shared_datadir, add_cli_arguments
):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "dummy_sani"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
            "remove_unsanitizable_products": True,
        }
    )
    add_cli_arguments(config_path)

    expansion_main()

    with open(shared_datadir / "dummy_sani_template_library.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 10


def test_preprocess_recommender(write_yaml, shared_datadir, add_cli_arguments):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "dummy"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
        }
    )
    add_cli_arguments(config_path)

    expansion_main()

    with open(shared_datadir / "dummy_template_library.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 10

    os.remove(shared_datadir / "dummy_training.csv")
    os.remove(shared_datadir / "dummy_testing.csv")
    os.remove(shared_datadir / "dummy_validation.csv")
    os.remove(shared_datadir / "dummy_unique_templates.hdf5")

    recommender_main()

    with open(shared_datadir / "dummy_training.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 6

    with open(shared_datadir / "dummy_testing.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2

    with open(shared_datadir / "dummy_validation.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2

    data = pd.read_hdf(shared_datadir / "dummy_unique_templates.hdf5", "table")
    assert len(data) == 2


def test_preprocess_filter(write_yaml, shared_datadir, add_cli_arguments):
    def duplicate_file(filename):
        with open(shared_datadir / filename, "r") as fileobj:
            lines = fileobj.read().splitlines()
        lines = lines + lines
        with open(shared_datadir / filename, "w") as fileobj:
            fileobj.write("\n".join(lines))

    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "make_false"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
            "library_headers": [
                "index",
                "reaction_hash",
                "reactants",
                "products",
                "retro_template",
                "template_hash",
                "template_code",
            ],
        }
    )

    add_cli_arguments(f"{config_path} strict")
    make_false_main()

    duplicate_file("make_false_template_library_false.csv")
    duplicate_file("make_false_template_library.csv")

    add_cli_arguments(config_path)
    filter_main()

    with open(shared_datadir / "make_false_training.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 6

    with open(shared_datadir / "make_false_testing.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2

    with open(shared_datadir / "make_false_validation.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2


def test_make_false_products(write_yaml, shared_datadir, add_cli_arguments):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "make_false"),
            "library_headers": [
                "index",
                "reaction_hash",
                "reactants",
                "products",
                "retro_template",
                "template_hash",
                "template_code",
            ],
        }
    )
    add_cli_arguments(f"{config_path} strict")

    make_false_main()

    with open(shared_datadir / "make_false_template_library_false.csv", "r") as fileobj:
        lines = fileobj.read().splitlines()
    assert len(lines) == 2


def test_download_public_data(tmpdir, mocker, add_cli_arguments):
    request_mock = mocker.patch("aizynthfinder.tools.download_public_data.requests.get")
    response_mock = request_mock.return_value
    filecontent = [b"abc", b"def"]
    response_mock.__enter__.return_value.iter_content.return_value = filecontent
    add_cli_arguments(str(tmpdir))

    download_main()

    filenames = glob.glob(str(tmpdir / "*.hdf5"))
    assert len(filenames) == 3
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
