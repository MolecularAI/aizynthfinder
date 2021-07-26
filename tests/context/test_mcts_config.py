import pytest

from aizynthfinder.context.config import Configuration
from aizynthfinder.aizynthfinder import AiZynthFinder


def test_load_empty_dict(default_config):
    config = Configuration.from_dict({})

    assert config == default_config


def test_load_empty_file(default_config, write_yaml):
    filename = write_yaml({})

    config = Configuration.from_file(filename)

    assert config == default_config


def test_load_from_file(write_yaml):
    filename = write_yaml(
        {
            "properties": {"cutoff_number": 300},
            "policy": {"properties": {"C": 1.9, "time_limit": 200}},
            "finder": {"properties": {"time_limit": 300}},
        }
    )

    config = Configuration.from_file(filename)

    assert config.cutoff_number == 300
    assert config.C == 1.9
    assert config.time_limit == 200


def test_load_from_dict_invalid_property(write_yaml):
    dict_ = {
        "properties": {"dummy": 300},
        "policy": {"properties": {"C": 1.9}},
        "finder": {"properties": {"time_limit": 300}},
    }
    with pytest.raises(AttributeError):
        Configuration.from_dict(dict_)


def test_get_properties(default_config):
    config = default_config

    props = config.properties

    # Just check a few properties
    assert props["C"] == 1.4
    assert props["time_limit"] == 120
    assert props["max_transforms"] == 6


def test_update_properties(default_config):
    config = default_config

    assert config.C != 2.0
    assert config.time_limit != 300
    assert config.max_transforms == 6

    config.properties = {"C": 2.0, "time_limit": 300, "max_transforms": None}

    assert config.C == 2.0
    assert config.time_limit == 300
    assert config.max_transforms == 6

    with pytest.raises(AttributeError):
        config.properties = {
            "C": 2.0,
            "time_limit": 300,
            "max_transforms": None,
            "dummy": 2,
        }


def test_load_stock(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml({"stock": {"files": {"test": stock_filename}}})

    config = Configuration.from_file(filename)

    assert config.stock.items == ["test"]


def test_load_policy(write_yaml, create_dummy_templates, mock_keras_model):
    templates_filename = create_dummy_templates(3)
    filename = write_yaml(
        {"policy": {"files": {"test": ["dummy", templates_filename]}}}
    )

    config = Configuration.from_file(filename)

    assert config.expansion_policy.items == ["test"]


def test_load_filter_policy(write_yaml, mock_keras_model):
    filename = write_yaml({"filter": {"files": {"test": "dummy"}}})

    config = Configuration.from_file(filename)

    assert config.filter_policy.items == ["test"]


def test_load_stop_criteria(write_yaml):
    filename = write_yaml(
        {"stock": {"stop_criteria": {"price": 100, "counts": {"C": 10}}}}
    )

    config = Configuration.from_file(filename)

    set_keys = [key for key, item in config.stock.stop_criteria.items() if item]
    assert set_keys == ["price", "counts"]
    assert list(config.stock.stop_criteria["counts"].keys()) == ["C"]
    assert config.stock.stop_criteria["counts"]["C"] == 10


def test_load_default_mongodb(write_yaml, mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.queries.get_mongo_client")
    filename = write_yaml({"stock": {"mongodb": {}}})

    config = Configuration.from_file(filename)

    mocked_client.assert_called_with("localhost")
    assert config.stock.items == ["mongodb_stock"]


def test_load_specific_mongodb(write_yaml, mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.queries.get_mongo_client")
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
    assert config.stock.items == ["mongodb_stock"]
    config.stock["mongodb_stock"].client.__getitem__.assert_called_with("mydatabase")
    config.stock["mongodb_stock"].database.__getitem__.assert_called_with(
        "mycollection"
    )


def test_load_external_stock(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.context.stock.stock.InMemoryInchiKeyQuery": {
                    "filename": stock_filename
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == ["InMemoryInchiKeyQuery"]


def test_load_external_stock_incorrect_module(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.context.stocks.stock.InMemoryInchiKeyQuery": {
                    "filename": stock_filename
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == []


def test_load_external_stock_incorrect_class(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "aizynthfinder.context.stock.InnMemoryInchiKeyQuery": {
                    "filename": stock_filename
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == []


def test_load_scorer_from_context_module(write_yaml):
    filename = write_yaml({"scorer": {"PriceSumScorer": None}})

    config = Configuration.from_file(filename)

    assert "sum of prices" in config.scorers.items


def test_load_scorer_from_module_spec(write_yaml):
    filename = write_yaml(
        {"scorer": {"aizynthfinder.context.scoring.PriceSumScorer": None}}
    )

    config = Configuration.from_file(filename)

    assert "sum of prices" in config.scorers.items


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
