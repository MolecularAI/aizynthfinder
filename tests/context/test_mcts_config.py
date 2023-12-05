import os
from unittest import mock

import pytest

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.context.config import Configuration
from aizynthfinder.context.stock import StockException


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
            "search": {
                "algorithm": "mcts",
                "algorithm_config": {
                    "C": 1.9,
                    "default_prior": 0.9,
                    "use_prior": False,
                },
                "max_transforms": 6,
                "time_limit": 200,
                "break_bonds_operator": "OR",
            },
        }
    )

    config = Configuration.from_file(filename)

    assert config.search.algorithm == "mcts"
    assert config.search.time_limit == 200

    assert config.search.algorithm_config == {
        "C": 1.9,
        "default_prior": 0.9,
        "use_prior": False,
        "prune_cycles_in_search": True,
        "search_reward": "state score",
        "immediate_instantiation": (),
        "mcts_grouping": None,
    }


def test_load_from_dict_invalid_property(write_yaml):
    dict_ = {
        "search": {"dummy": 300, "algorithm": "mcts"},
    }
    with pytest.raises(AttributeError, match="Could not find attribute to set: dummy"):
        Configuration.from_dict(dict_)


def test_update_search(default_config):
    config = default_config

    assert config.search.algorithm_config["C"] != 2.0
    assert config.search.time_limit != 300
    assert config.search.max_transforms == 6

    config = config.from_dict(
        {
            "search": {
                "algorithm_config": {"C": 2.0},
                "time_limit": 300,
                "max_transforms": None,
            }
        }
    )
    assert config.search.algorithm_config["C"] == 2.0
    assert config.search.time_limit == 300
    assert config.search.max_transforms == 6

    with pytest.raises(AttributeError):
        config = config.from_dict(
            {
                "search": {
                    "time_limit": 300,
                    "max_transforms": None,
                    "dummy": 2,
                }
            }
        )


def test_load_stock(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {"stock": {"buyables": {"type": "inchiset", "path": stock_filename}}}
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == ["test"]


def test_load_policy(write_yaml, create_dummy_templates, mock_onnx_model):
    templates_filename = create_dummy_templates(3)
    filename = write_yaml(
        {
            "expansion": {
                "uspto": {
                    "type": "template-based",
                    "model": "dummy.onnx",
                    "template": templates_filename,
                    "cutoff_number": 75,
                },
                "full": ["dummy.onnx", templates_filename],
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.expansion_policy.items == ["full", "uspto"]


def test_load_filter_policy(write_yaml, mock_onnx_model):
    filename = write_yaml(
        {
            "filter": {
                "uspto": {
                    "type": "quick-filter",
                    "model": "dummy.onnx",
                    "exclude_from_policy": ["rc"],
                },
                "full": "dummy.onnx",
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.filter_policy.items == ["full", "uspto"]


def test_load_stock(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "buyables": {"type": "inchiset", "path": stock_filename},
                "emolecules": stock_filename,
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == ["buyables", "emolecules"]


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
    filename = write_yaml({"stock": {"mongodb_stock": {"type": "mongodb"}}})

    config = Configuration.from_file(filename)

    mocked_client.assert_called_with("localhost")
    assert config.stock.items == ["mongodb_stock"]


def test_load_specific_mongodb(write_yaml, mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.queries.get_mongo_client")
    filename = write_yaml(
        {
            "stock": {
                "mongodb_stock": {
                    "type": "mongodb",
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
                "inchi": {
                    "type": "aizynthfinder.context.stock.stock.InMemoryInchiKeyQuery",
                    "path": stock_filename,
                }
            }
        }
    )

    config = Configuration.from_file(filename)

    assert config.stock.items == ["inchi"]


def test_load_external_stock_incorrect_module(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "inchi": {
                    "type": "aizynthfinder.context.stocks.stock.InMemoryInchiKeyQuery",
                    "path": stock_filename,
                }
            }
        }
    )

    with pytest.raises(StockException):
        config = Configuration.from_file(filename)


def test_load_external_stock_incorrect_class(write_yaml, create_dummy_stock1):
    stock_filename = create_dummy_stock1("hdf5")
    filename = write_yaml(
        {
            "stock": {
                "inchi": {
                    "type": "aizynthfinder.context.stock.InnMemoryInchiKeyQuery",
                    "path": stock_filename,
                }
            }
        }
    )

    with pytest.raises(StockException):
        config = Configuration.from_file(filename)


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


@mock.patch.dict(
    os.environ,
    {"ITERATION_LIMIT": "300", "C": "1.9"},
)
def test_load_yaml_with_environ(write_yaml):
    filename = write_yaml(
        {
            "search": {
                "algorithm_config": {"C": "${C}"},
                "iteration_limit": "${ITERATION_LIMIT}",
            },
        }
    )
    config = Configuration.from_file(filename)

    assert config.search.iteration_limit == 300
    assert config.search.algorithm_config["C"] == 1.9


@mock.patch.dict(os.environ, {"ITERATION_LIMIT": "300"})
def test_load_yaml_with_environ_raises_error(write_yaml):
    filename = write_yaml(
        {
            "search": {
                "algorithm_config": {"C": 1.9},
                "iteration_limit": "${ITERATION}",
            },
        }
    )
    with pytest.raises(ValueError, match="'ITERATION' not in environment variables"):
        Configuration.from_file(filename)


def test_load_algorithm_config(write_yaml):
    filename = write_yaml({"search": {"algorithm_config": {"C": 1.9}}})

    config = Configuration.from_file(filename)

    expected_keys = [
        "C",
        "default_prior",
        "use_prior",
        "prune_cycles_in_search",
        "search_reward",
        "immediate_instantiation",
        "mcts_grouping",
    ]
    for key in expected_keys:
        assert key in config.search.algorithm_config, f"{key} not in config"
    assert config.search.algorithm_config["C"] == 1.9


def test_load_algorithm_config_failure(write_yaml):
    filename = write_yaml({"search": {"algorithm_config": 5.5}})

    with pytest.raises(ValueError, match="algorithm_config"):
        Configuration.from_file(filename)


def test_init_search_yaml(write_yaml, create_dummy_templates, mock_onnx_model):
    templates_filename = create_dummy_templates(3)
    filename = write_yaml(
        {
            "expansion": {
                "policy1": {
                    "type": "template-based",
                    "model": "dummy.onnx",
                    "template": templates_filename,
                    "cutoff_number": 300,
                }
            },
            "search": {"time_limit": 300},
        }
    )

    finder = AiZynthFinder(filename)

    assert finder.config.expansion_policy["policy1"].cutoff_number == 300
    assert finder.config.search.time_limit == 300


def test_init_search_dict(create_dummy_templates, mock_onnx_model):
    templates_filename = create_dummy_templates(3)
    dict_ = {
        "expansion": {
            "policy1": {
                "type": "template-based",
                "model": "dummy.onnx",
                "template": templates_filename,
                "cutoff_number": 300,
            }
        },
        "search": {"time_limit": 300},
    }

    finder = AiZynthFinder(configdict=dict_)

    assert finder.config.expansion_policy["policy1"].cutoff_number == 300
    assert finder.config.search.time_limit == 300


def test_init_search_yaml_dict(write_yaml, create_dummy_templates, mock_onnx_model):
    templates_filename = create_dummy_templates(3)
    filename = write_yaml(
        {
            "expansion": {
                "policy1": {
                    "type": "template-based",
                    "model": "dummy.onnx",
                    "template": templates_filename,
                    "cutoff_number": 300,
                }
            },
            "search": {"time_limit": 300},
        }
    )
    dict_ = {
        "expansion": {
            "policy1": {
                "type": "template-based",
                "model": "dummy.onnx",
                "template": templates_filename,
                "cutoff_number": 100,
            }
        },
    }

    finder = AiZynthFinder(filename, configdict=dict_)

    assert finder.config.expansion_policy["policy1"].cutoff_number == 300


def test_init_search_none(default_config):
    finder = AiZynthFinder()

    assert finder.config == default_config
