import pandas as pd

from aizynthfinder.interfaces import AiZynthApp
from aizynthfinder.interfaces.aizynthapp import main as app_main
from aizynthfinder.interfaces.aizynthcli import main as cli_main
from aizynthfinder.tools.make_stock import main as make_stock_main
from aizynthfinder.training.preprocess_rollout import main as rollout_main
from aizynthfinder.training.utils import Config


def test_create_gui_app(mocker):
    display_patch = mocker.patch("aizynthfinder.interfaces.aizynthapp.display")
    AiZynthApp(configfile=None, setup=False)

    display_patch.assert_not_called()

    AiZynthApp(configfile=None)

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

    stock.load_stock(filename, "stock1")
    stock.select_stocks(["stock1"])
    assert len(stock) == 3


def test_preprocess_rollout(write_yaml, shared_datadir, add_cli_arguments):
    config_path = write_yaml(
        {
            "file_prefix": str(shared_datadir / "dummy"),
            "split_size": {"training": 0.6, "testing": 0.2, "validation": 0.2},
        }
    )
    add_cli_arguments(config_path)

    rollout_main()

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
