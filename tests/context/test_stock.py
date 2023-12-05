import pytest
import sys

import pandas as pd

from aizynthfinder.chem import Molecule
from aizynthfinder.context.stock import (
    StockException,
)
from aizynthfinder.context.stock.queries import HAS_MOLBLOOM
from aizynthfinder.tools.make_stock import (
    extract_plain_smiles,
    extract_smiles_from_module,
    make_hdf5_stock,
    make_mongo_stock,
    make_molbloom,
    make_molbloom_inchi,
)


def test_load_stock(default_config, setup_stock_with_query):
    stock_query = setup_stock_with_query()

    default_config.stock.load(stock_query, "stock1")

    assert "stock1" in default_config.stock.items

    with pytest.raises(StockException):
        default_config.stock.load("stock", "stock1")


def test_load_csv_stock(default_config, create_dummy_stock1, setup_stock_with_query):
    stock_query = setup_stock_with_query(create_dummy_stock1("csv"))

    default_config.stock.load(stock_query, "stock1")

    assert len(default_config.stock["stock1"]) == 2


def test_load_txt_stock(default_config, create_dummy_stock1, setup_stock_with_query):
    stock_query = setup_stock_with_query(create_dummy_stock1("txt"))

    default_config.stock.load(stock_query, "stock1")

    assert len(default_config.stock["stock1"]) == 2


def test_load_from_config(default_config, create_dummy_stock1):
    filename = create_dummy_stock1("csv")
    stock = default_config.stock

    stock.load_from_config(
        **{
            "stock1": {
                "type": "InMemoryInchiKeyQuery",
                "path": filename,
                "price_col": "price",
            },
            "stock2": filename,
            "stock3": {
                "type": "aizynthfinder.context.stock.queries.InMemoryInchiKeyQuery",
                "path": filename,
            },
            "stock4": {
                "type": "inchiset",
                "path": filename,
            },
        }
    )

    assert "stock1" in stock.items
    assert "stock2" in stock.items
    assert "stock3" in stock.items
    assert "stock4" in stock.items


def test_availability_string(
    default_config, setup_stock_with_query, create_dummy_stock1, create_dummy_stock2
):
    stock = default_config.stock
    stock.load(setup_stock_with_query(create_dummy_stock1("hdf5")), "stock1")
    stock.load(setup_stock_with_query(create_dummy_stock2), "stock2")

    stock.select(["stock1", "stock2"])

    ethanol = Molecule(smiles="CCO")
    assert stock.availability_string(ethanol) == "Not in stock"

    benzene = Molecule(smiles="c1ccccc1")
    assert stock.availability_string(benzene) == "stock1,stock2"

    toluene = Molecule(smiles="Cc1ccccc1")
    assert stock.availability_string(toluene) == "stock1"


def test_mol_in_stock(setup_stock_with_query):
    stock = setup_stock_with_query()

    ethanol = Molecule(smiles="CCO")
    assert ethanol not in stock

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    toluene = Molecule(smiles="Cc1ccccc1")
    assert toluene in stock


def test_load_csv_stock_with_price(default_config, create_dummy_stock1):
    filename = create_dummy_stock1("csv")

    default_config.stock.load_from_config(
        **{
            "stock1": {
                "type": "InMemoryInchiKeyQuery",
                "path": filename,
                "price_col": "price",
            }
        }
    )
    default_config.stock.select_all()

    assert len(default_config.stock) == 2

    benzene = Molecule(smiles="c1ccccc1")
    assert default_config.stock.price(benzene) == 5.0

    toluene = Molecule(smiles="Cc1ccccc1")
    assert default_config.stock.price(toluene) == 10.0


def test_load_csv_stock_with_price_duplicates(tmpdir, default_config):
    filename = str(tmpdir / "stock1.hdf5")
    pd.DataFrame(
        {
            "inchi_key": ["UHOVQNZJYSORNB-UHFFFAOYSA-N", "UHOVQNZJYSORNB-UHFFFAOYSA-N"],
            "price": [5.0, 10.0],
        }
    ).to_hdf(filename, "table")

    with pytest.raises(StockException, match="unique"):
        default_config.stock.load_from_config(
            **{
                "stock1": {
                    "type": "InMemoryInchiKeyQuery",
                    "path": filename,
                    "price_col": "price",
                }
            }
        )


def test_load_csv_stock_with_null_price(tmpdir, default_config):
    filename = str(tmpdir / "stock1.hdf5")
    pd.DataFrame(
        {
            "inchi_key": ["UHOVQNZJYSORNB-UHFFFAOYSA-N", "YXFVVABEGXRONW-UHFFFAOYSA-N"],
            "price": [5.0, None],
        }
    ).to_hdf(filename, "table")

    with pytest.raises(StockException, match="impute"):
        default_config.stock.load_from_config(
            **{
                "stock1": {
                    "type": "InMemoryInchiKeyQuery",
                    "path": filename,
                    "price_col": "price",
                }
            }
        )


def test_load_csv_stock_with_negative_price(tmpdir, default_config):
    filename = str(tmpdir / "stock1.hdf5")
    pd.DataFrame(
        {
            "inchi_key": ["UHOVQNZJYSORNB-UHFFFAOYSA-N", "YXFVVABEGXRONW-UHFFFAOYSA-N"],
            "price": [5.0, -1],
        }
    ).to_hdf(filename, "table")

    with pytest.raises(StockException, match="non-negative"):
        default_config.stock.load_from_config(
            **{
                "stock1": {
                    "type": "InMemoryInchiKeyQuery",
                    "path": filename,
                    "price_col": "price",
                }
            }
        )


def test_exclude(default_config, setup_stock_with_query):
    stock_query = setup_stock_with_query()
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    stock.exclude(benzene)
    assert benzene not in stock


def test_exclude_many(default_config, setup_stock_with_query):
    stock_query = setup_stock_with_query()
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    stock.exclude(benzene)
    assert benzene not in stock

    toluene = Molecule(smiles="Cc1ccccc1")
    assert toluene in stock

    stock.exclude(toluene)
    assert benzene not in stock
    assert toluene not in stock

    stock.reset_exclusion_list()
    stock.exclude(toluene)
    assert benzene in stock
    assert toluene not in stock


def test_price_no_price(default_config, setup_stock_with_query):
    stock_query = setup_stock_with_query()
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_price_with_price(default_config, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    price = {mol: 14}
    stock_query = make_stock_query([mol], price=price)
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    assert stock.price(mol) == 14


def test_price_with_price_raises(default_config, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock_query = make_stock_query([mol])
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_amount_no_amount(default_config, setup_stock_with_query):
    stock_query = setup_stock_with_query()
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_amount_with_amount(default_config, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    amount = {mol: 14}
    stock_query = make_stock_query([mol], amount=amount)
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    assert stock.amount(mol) == 14


def test_amount_with_amount_raises(default_config, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock_query = make_stock_query([mol])
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_counts_filter(default_config, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    stock_query = make_stock_query([mol1, mol2])
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"size": {"C": 2}})

    assert mol1 not in stock

    stock.set_stop_criteria({"size": {"C": 10}})

    assert mol1 in stock

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})

    assert mol1 in stock
    assert mol2 not in stock


def test_no_entries_filter(default_config, make_stock_query):
    stock_query = make_stock_query([])
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})
    assert Molecule(smiles="c1ccccc1") not in stock


def test_amounts_filter(default_config, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    stock_query = make_stock_query([mol1, mol2, mol3], amount={mol1: 10, mol2: 5})
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"amount": 10})

    assert mol1 in stock
    assert mol2 not in stock
    assert mol3 in stock


def test_price_filter(default_config, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    stock_query = make_stock_query([mol1, mol2, mol3], price={mol1: 10, mol2: 5})
    stock = default_config.stock

    stock.load(stock_query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"price": 5})

    assert mol1 not in stock
    assert mol2 in stock
    assert mol3 in stock


def test_mongodb_load_default(mocked_mongo_db_query):
    mocked_client, query = mocked_mongo_db_query()

    mocked_client.assert_called_with("localhost")
    query.client.__getitem__.assert_called_with("stock_db")
    query.database.__getitem__.assert_called_with("molecules")


def test_mongodb_load_specifics(mocked_mongo_db_query):
    mocked_client, query = mocked_mongo_db_query(
        host="127.0.0", database="my_database", collection="my_collection"
    )

    mocked_client.assert_called_with("127.0.0")
    query.client.__getitem__.assert_called_with("my_database")
    query.database.__getitem__.assert_called_with("my_collection")


def test_mongodb_load_host_env(mocked_mongo_db_query, monkeypatch):
    monkeypatch.setenv("MONGODB_HOST", "myhost")

    mocked_client, query = mocked_mongo_db_query()

    mocked_client.assert_called_with("myhost")
    query.client.__getitem__.assert_called_with("stock_db")
    query.database.__getitem__.assert_called_with("molecules")


def test_mongodb_contains(mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.count_documents.side_effect = [1, 0]

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in query

    toluene = Molecule(smiles="Cc1ccccc1")
    assert toluene not in query


def test_mongodb_availability(mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.find.return_value = [{"source": "source1"}, {"source": "source2"}]

    benzene = Molecule(smiles="c1ccccc1")
    assert query.availability_string(benzene) == "source1,source2"


def test_mongodb_integration(default_config, mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.count_documents.return_value = 1
    stock = default_config.stock

    stock.load(query, "stock1")
    stock.select(["stock1"])

    ethanol = Molecule(smiles="CCO")
    assert ethanol in stock


def test_extract_smiles_from_plain_file(create_dummy_smiles_source):
    filename = create_dummy_smiles_source("txt")

    smiles = list(extract_plain_smiles([filename, filename]))

    assert len(smiles) == 8


def test_extract_smiles_from_module(create_dummy_smiles_source, shared_datadir):
    filename = create_dummy_smiles_source("csv")
    module_name = "custom_loader"
    sys.path.append(str(shared_datadir))

    smiles = list(extract_smiles_from_module([module_name, filename, filename]))

    assert len(smiles) == 8

    module_name = "custom_loader2"
    smiles = list(extract_smiles_from_module([module_name]))

    assert len(smiles) == 4

    sys.path.pop()


def test_make_hdf5_stock(default_config, tmpdir):
    filename = str(tmpdir / "temp.hdf5")
    inchi_keys = ("key1", "key2", "key1")

    make_hdf5_stock(inchi_keys, filename)
    stock = default_config.stock
    stock.load_from_config(
        **{
            "stock1": {
                "type": "inchiset",
                "path": filename,
            }
        }
    )
    stock.select(["stock1"])

    assert len(stock) == 2


def test_make_mongodb_stock(mocked_mongo_db_query):
    inchi_keys = ("key1", "key2", "key1")
    _, query = mocked_mongo_db_query()

    make_mongo_stock(inchi_keys, "temp")

    query.molecules.insert_many.assert_called_once()
    # [0][0] is to get the first ordered argument, i.e. the generator of items to be inserted
    assert len(list(query.molecules.insert_many.call_args[0][0])) == 2


@pytest.mark.xfail(condition=not HAS_MOLBLOOM, reason="molbloom package not installed")
def test_load_bloom_filter(default_config, shared_datadir):
    stock = default_config.stock
    filename = str(shared_datadir / "simple_filter.bloom")

    stock.load_from_config(**{"bloom": filename})

    assert Molecule(smiles="c1ccccc1") in stock["bloom"]


@pytest.mark.xfail(condition=not HAS_MOLBLOOM, reason="molbloom package not installed")
def test_make_bloom_filter(default_config, tmpdir):
    stock = default_config.stock
    filename = str(tmpdir / "temp.bloom")
    smiles = ("CC", "CCO")

    make_molbloom(smiles, filename, 10, 3)

    stock.load_from_config(molbloom=filename)

    assert Molecule(smiles="CC") in stock["molbloom"]


@pytest.mark.xfail(condition=not HAS_MOLBLOOM, reason="molbloom package not installed")
def test_make_bloom_inchi_filter(default_config, tmpdir):
    stock = default_config.stock
    filename = str(tmpdir / "temp2.bloom")
    inchi_keys = ("label1", "label2")

    make_molbloom_inchi(inchi_keys, filename, 10, 3)

    stock.load_from_config(molbloom=filename)

    assert "molbloom" in stock.items
