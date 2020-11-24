from aizynthfinder.mcts.config import Configuration
from aizynthfinder.aizynthfinder import AiZynthFinder


def test_load_empty_dict(default_config):
    config = Configuration.from_dict({})

    assert config == default_config


def test_load_empty_file(default_config, write_yaml):
    filename = write_yaml({})

    config = Configuration.from_file(filename)

    assert config == default_config


def test_update_properties(write_yaml):
    filename = write_yaml(
        {
            "properties": {"cutoff_number": 300},
            "policy": {"properties": {"C": 1.9}},
            "finder": {"properties": {"time_limit": 300}},
        }
    )

    config = Configuration.from_file(filename)

    assert config.cutoff_number == 300
    assert config.C == 1.9
    assert config.time_limit == 300


def test_load_stock(write_yaml, shared_datadir):
    stock_filename = str(shared_datadir / "stock1.h5")
    filename = write_yaml({"stock": {"files": {"test": stock_filename}}})

    config = Configuration.from_file(filename)

    assert config.stock.available_stocks() == ["test"]


def test_load_policy(write_yaml, shared_datadir, mock_policy_model):
    templates_filename = str(shared_datadir / "templates.hdf5")
    filename = write_yaml(
        {"policy": {"files": {"test": ["dummy", templates_filename]}}}
    )

    config = Configuration.from_file(filename)

    assert config.policy.available_policies() == ("test",)


def test_load_default_mongodb(write_yaml, mocker):
    mocked_client = mocker.patch("aizynthfinder.mcts.stock.MongoClient")
    filename = write_yaml({"stock": {"mongodb": {}}})

    config = Configuration.from_file(filename)

    mocked_client.assert_called_with(None)
    assert config.stock.available_stocks() == ["mongodb_stock"]


def test_load_specific_mongodb(write_yaml, mocker):
    mocked_client = mocker.patch("aizynthfinder.mcts.stock.MongoClient")
    filename = write_yaml(
        {
            "stock": {
                "mongodb": {
                    "host": "myhost",
                    "database": "mydatabase",
                    "collection": "mycollection",
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    mocked_client.assert_called_with("myhost")
    assert config.stock.available_stocks() == ["mongodb_stock"]
    config.stock["mongodb_stock"].client.__getitem__.assert_called_with("mydatabase")
    config.stock["mongodb_stock"].database.__getitem__.assert_called_with(
        "mycollection"
    )


def test_load_external_stock(write_yaml, shared_datadir):
    stock_filename = str(shared_datadir / "stock1.h5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.mcts.stock.InMemoryInchiKeyQuery": {
                    "filename": stock_filename
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.available_stocks() == ["InMemoryInchiKeyQuery"]


def test_load_external_stock_incorrect_module(write_yaml, shared_datadir):
    stock_filename = str(shared_datadir / "stock1.h5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.mcts.InMemoryInchiKeyQuery": {"filename": stock_filename}
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.available_stocks() == []


def test_load_external_stock_incorrect_class(write_yaml, shared_datadir):
    stock_filename = str(shared_datadir / "stock1.h5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.mcts.stock.InnMemoryInchiKeyQuery": {
                    "filename": stock_filename
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.available_stocks() == []


def test_init_search_yaml(write_yaml):
    filename = write_yaml(
        {
            "properties": {"cutoff_number": 300},
            "policy": {"properties": {"C": 1.9}},
            "finder": {"properties": {"time_limit": 300}},
        }
    )

    finder = AiZynthFinder(filename)

    assert finder.config.cutoff_number == 300
    assert finder.config.C == 1.9
    assert finder.config.time_limit == 300


def test_init_search_dict():
    dict_ = {
        "properties": {"cutoff_number": 300},
        "policy": {"properties": {"C": 1.9}},
        "finder": {"properties": {"time_limit": 300}},
    }

    finder = AiZynthFinder(configdict=dict_)

    assert finder.config.cutoff_number == 300
    assert finder.config.C == 1.9
    assert finder.config.time_limit == 300


def test_init_search_yaml_dict(write_yaml):
    filename = write_yaml({"properties": {"cutoff_number": 300}})
    dict_ = {
        "properties": {"cutoff_number": 100},
    }

    finder = AiZynthFinder(filename, configdict=dict_)

    assert finder.config.cutoff_number == 300


def test_init_search_none(default_config):

    finder = AiZynthFinder()

    assert finder.config == default_config
