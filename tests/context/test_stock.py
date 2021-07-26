import pytest
import sys

from aizynthfinder.chem import Molecule
from aizynthfinder.context.stock import (
    StockException,
)
from aizynthfinder.tools.make_stock import (
    extract_plain_smiles,
    extract_smiles_from_module,
    make_hdf5_stock,
    make_mongo_stock,
)


def test_load_stock(default_config, create_dummy_stock1):
    filename = create_dummy_stock1("hdf5")

    default_config.stock.load(filename, "stock1")

    assert len(default_config.stock["stock1"]) == 2


def test_load_csv_stock(default_config, create_dummy_stock1):
    filename = create_dummy_stock1("csv")

    default_config.stock.load(filename, "stock1")

    assert len(default_config.stock["stock1"]) == 2


def test_load_txt_stock(default_config, create_dummy_stock1):
    filename = create_dummy_stock1("txt")

    default_config.stock.load(filename, "stock1")

    assert len(default_config.stock["stock1"]) == 2


def test_availability_string(default_config, create_dummy_stock1, create_dummy_stock2):
    stock = default_config.stock
    stock.load(create_dummy_stock1("hdf5"), "stock1")
    stock.load(create_dummy_stock2, "stock2")

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


def test_exclude(setup_stock_with_query):
    stock = setup_stock_with_query()

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    stock.exclude(benzene)
    assert benzene not in stock


def test_exclude_many(setup_stock_with_query):
    stock = setup_stock_with_query()

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


def test_price_no_price(setup_stock_with_query):
    stock = setup_stock_with_query()

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_price_with_price(setup_stock_with_query, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    price = {mol: 14}
    stock = setup_stock_with_query(make_stock_query([mol], price=price))

    assert stock.price(mol) == 14


def test_price_with_price_raises(setup_stock_with_query, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock = setup_stock_with_query(make_stock_query([mol]))

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_amount_no_amount(setup_stock_with_query):
    stock = setup_stock_with_query()

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_amount_with_amount(setup_stock_with_query, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    amount = {mol: 14}
    stock = setup_stock_with_query(make_stock_query([mol], amount=amount))

    assert stock.amount(mol) == 14


def test_amount_with_amount_raises(setup_stock_with_query, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock = setup_stock_with_query(make_stock_query([mol]))

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_counts_filter(setup_stock_with_query, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    stock = setup_stock_with_query(make_stock_query([mol1, mol2]))

    stock.set_stop_criteria({"size": {"C": 2}})

    assert mol1 not in stock

    stock.set_stop_criteria({"size": {"C": 10}})

    assert mol1 in stock

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})

    assert mol1 in stock
    assert mol2 not in stock


def test_no_entries_filter(setup_stock_with_query, make_stock_query):
    stock = setup_stock_with_query(make_stock_query([]))

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})
    assert Molecule(smiles="c1ccccc1") not in stock


def test_amounts_filter(setup_stock_with_query, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    query = make_stock_query([mol1, mol2, mol3], amount={mol1: 10, mol2: 5})
    stock = setup_stock_with_query(query)

    stock.set_stop_criteria({"amount": 10})

    assert mol1 in stock
    assert mol2 not in stock
    assert mol3 in stock


def test_price_filter(setup_stock_with_query, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    query = make_stock_query([mol1, mol2, mol3], price={mol1: 10, mol2: 5})
    stock = setup_stock_with_query(query)

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


def test_mongodb_integration(setup_stock_with_query, mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.count_documents.return_value = 1
    stock = setup_stock_with_query(query)

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
    stock.load(filename, "stock1")
    stock.select(["stock1"])

    assert len(stock) == 2


def test_make_mongodb_stock(mocked_mongo_db_query):
    inchi_keys = ("key1", "key2", "key1")
    _, query = mocked_mongo_db_query()

    make_mongo_stock(inchi_keys, "temp")

    query.molecules.insert_many.assert_called_once()
    # [0][0] is to get the first ordered argument, i.e. the generator of items to be inserted
    assert len(list(query.molecules.insert_many.call_args[0][0])) == 2
