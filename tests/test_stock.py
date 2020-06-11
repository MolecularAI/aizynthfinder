import pytest
import sys

from aizynthfinder.chem import Molecule
from aizynthfinder.mcts.stock import StockException, MongoDbInchiKeyQuery
from aizynthfinder.tools.make_stock import (
    extract_plain_smiles,
    extract_smiles_from_module,
    make_hdf5_stock,
    make_mongo_stock,
)


def test_load_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")

    stock.load_stock(filename, "stock1")

    assert len(stock) == 0


def test_load_csv_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.csv")

    stock.load_stock(filename, "stock1")

    assert len(stock) == 0


def test_load_txt_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.txt")

    stock.load_stock(filename, "stock1")

    assert len(stock) == 0


def test_select_one_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")

    stock.select_stocks(["stock1"])

    assert stock["stock1"].stock_inchikeys == {
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "YXFVVABEGXRONW-UHFFFAOYSA-N",
    }


def test_select_one_stock_as_str(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")

    stock.select_stocks("stock1")

    assert stock["stock1"].stock_inchikeys == {
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "YXFVVABEGXRONW-UHFFFAOYSA-N",
    }


def test_select_one_stock_load_two(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    filename = str(shared_datadir / "stock2.h5")
    stock.load_stock(filename, "stock2")

    stock.select_stocks(["stock1"])

    assert stock["stock1"].stock_inchikeys == {
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "YXFVVABEGXRONW-UHFFFAOYSA-N",
    }

    assert stock["stock2"].stock_inchikeys == {
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
        "ISWSIDIOOBJBQZ-UHFFFAOYSA-N",
    }

    assert len(stock) == 2


def test_select_two_stocks(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    filename = str(shared_datadir / "stock2.h5")
    stock.load_stock(filename, "stock2")

    stock.select_stocks(["stock1", "stock2"])

    assert len(stock) == 4


def test_select_invalid_key(stock, shared_datadir, mocker):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")

    with pytest.raises(StockException):
        stock.select_stocks(["stock3", "stock1"])


def test_select_and_add_remove(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    filename = str(shared_datadir / "stock2.h5")
    stock.load_stock(filename, "stock2")

    stock.select_stocks(["stock1"])
    stock.add_stock("stock2")

    assert len(stock) == 4

    stock.remove_stock("stock2")

    assert len(stock) == 2


def test_availability_string(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    filename = str(shared_datadir / "stock2.h5")
    stock.load_stock(filename, "stock2")

    stock.select_stocks(["stock1", "stock2"])

    ethanol = Molecule(smiles="CCO")
    assert stock.availability_string(ethanol) == "Not in stock"

    benzene = Molecule(smiles="c1ccccc1")
    assert stock.availability_string(benzene) == "stock1,stock2"

    toluene = Molecule(smiles="Cc1ccccc1")
    assert stock.availability_string(toluene) == "stock1"


def test_mol_in_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    stock.select_stocks(["stock1"])

    ethanol = Molecule(smiles="CCO")
    assert ethanol not in stock

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    toluene = Molecule(smiles="Cc1ccccc1")
    assert toluene in stock


def test_exclude(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    stock.select_stocks(["stock1"])

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    stock.exclude(benzene)
    assert benzene not in stock


def test_exclude_many(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load_stock(filename, "stock1")
    stock.select_stocks(["stock1"])

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


@pytest.fixture
def mocked_mongo_db_query(mocker):
    mocked_client = mocker.patch("aizynthfinder.mcts.stock.MongoClient")

    def wrapper(**kwargs):
        return mocked_client, MongoDbInchiKeyQuery(**kwargs)

    return wrapper


def test_mongodb_load_default(mocked_mongo_db_query):
    mocked_client, query = mocked_mongo_db_query()

    mocked_client.assert_called_with(None)
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


def test_mongodb_integration(stock, mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.count_documents.return_value = 1
    stock.load_stock(query, "stock1")
    stock.select_stocks(["stock1"])

    ethanol = Molecule(smiles="CCO")
    assert ethanol in stock


def test_extract_smiles_from_plain_file(shared_datadir):
    filename = str(shared_datadir / "smiles_plain.txt")

    smiles = list(extract_plain_smiles([filename, filename]))

    assert len(smiles) == 8


def test_extract_smiles_from_module(shared_datadir):
    filename = str(shared_datadir / "smiles_csv.csv")
    module_name = "custom_loader"
    sys.path.append(str(shared_datadir))

    smiles = list(extract_smiles_from_module([module_name, filename, filename]))

    assert len(smiles) == 8

    module_name = "custom_loader2"
    smiles = list(extract_smiles_from_module([module_name]))

    assert len(smiles) == 4

    sys.path.pop()


def test_make_hdf5_stock(stock, tmpdir):
    filename = str(tmpdir / "temp.hdf5")
    inchi_keys = ("key1", "key2", "key1")

    make_hdf5_stock(inchi_keys, filename)
    stock.load_stock(filename, "stock1")
    stock.select_stocks(["stock1"])

    assert len(stock) == 2


def test_make_mongodb_stock(mocked_mongo_db_query):
    inchi_keys = ("key1", "key2", "key1")
    _, query = mocked_mongo_db_query()

    make_mongo_stock(inchi_keys, "temp")

    query.molecules.insert_many.assert_called_once()
    # [0][0] is to get the first ordered argument, i.e. the generator of items to be inserted
    assert len(list(query.molecules.insert_many.call_args[0][0])) == 2
