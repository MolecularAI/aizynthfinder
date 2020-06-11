from aizynthfinder.mcts.config import Configuration


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
