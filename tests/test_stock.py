import pytest
import sys

from aizynthfinder.chem import Molecule
from aizynthfinder.context.stock import (
    MongoDbInchiKeyQuery,
    StockQueryMixin,
    StockException,
)
from aizynthfinder.tools.make_stock import (
    extract_plain_smiles,
    extract_smiles_from_module,
    make_hdf5_stock,
    make_mongo_stock,
)


@pytest.fixture
def make_stock_query():
    class StockQuery(StockQueryMixin):
        def __init__(self, mols, price, amount):
            self.mols = mols
            self._price = price
            self._amount = amount

        def __contains__(self, mol):
            return mol in self.mols

        def __len__(self):
            return len(self.mols)

        def amount(self, mol):
            if mol in self._amount:
                return self._amount[mol]
            raise StockException()

        def price(self, mol):
            if mol in self._price:
                return self._price[mol]
            raise StockException()

    def wrapper(mols, amount=None, price=None):
        return StockQuery(mols, price or {}, amount or {})

    return wrapper


def test_load_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")

    stock.load(filename, "stock1")

    assert len(stock) == 0


def test_load_csv_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.csv")

    stock.load(filename, "stock1")

    assert len(stock) == 0


def test_load_txt_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.txt")

    stock.load(filename, "stock1")

    assert len(stock) == 0


def test_availability_string(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
    filename = str(shared_datadir / "stock2.h5")
    stock.load(filename, "stock2")

    stock.select(["stock1", "stock2"])

    ethanol = Molecule(smiles="CCO")
    assert stock.availability_string(ethanol) == "Not in stock"

    benzene = Molecule(smiles="c1ccccc1")
    assert stock.availability_string(benzene) == "stock1,stock2"

    toluene = Molecule(smiles="Cc1ccccc1")
    assert stock.availability_string(toluene) == "stock1"


def test_mol_in_stock(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
    stock.select(["stock1"])

    ethanol = Molecule(smiles="CCO")
    assert ethanol not in stock

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    toluene = Molecule(smiles="Cc1ccccc1")
    assert toluene in stock


def test_exclude(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
    stock.select(["stock1"])

    benzene = Molecule(smiles="c1ccccc1")
    assert benzene in stock

    stock.exclude(benzene)
    assert benzene not in stock


def test_exclude_many(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
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


def test_price_no_price(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_price_with_price(stock, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    price = {mol: 14}
    stock.load(make_stock_query([mol], price=price), "stock1")
    stock.select(["stock1"])

    assert stock.price(mol) == 14


def test_price_with_price_raises(stock, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock.load(make_stock_query([mol]), "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.price(Molecule(smiles="c1ccccc1"))


def test_amount_no_amount(stock, shared_datadir):
    filename = str(shared_datadir / "stock1.h5")
    stock.load(filename, "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_amount_with_amount(stock, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    amount = {mol: 14}
    stock.load(make_stock_query([mol], amount=amount), "stock1")
    stock.select(["stock1"])

    assert stock.amount(mol) == 14


def test_amount_with_amount_raises(stock, make_stock_query):
    mol = Molecule(smiles="c1ccccc1")
    stock.load(make_stock_query([mol]), "stock1")
    stock.select(["stock1"])

    with pytest.raises(StockException):
        stock.amount(Molecule(smiles="c1ccccc1"))


def test_counts_filter(stock, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    stock.load(make_stock_query([mol1, mol2]), "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"size": {"C": 2}})

    assert mol1 not in stock

    stock.set_stop_criteria({"size": {"C": 10}})

    assert mol1 in stock

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})

    assert mol1 in stock
    assert mol2 not in stock


def test_no_entries_filter(stock, make_stock_query):
    stock.load(make_stock_query([]), "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"size": {"C": 10, "O": 0}})

    assert Molecule(smiles="c1ccccc1") not in stock


def test_amounts_filter(stock, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    query = make_stock_query([mol1, mol2, mol3], amount={mol1: 10, mol2: 5})
    stock.load(query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"amount": 10})

    assert mol1 in stock
    assert mol2 not in stock
    assert mol3 in stock


def test_price_filter(stock, make_stock_query):
    mol1 = Molecule(smiles="c1ccccc1")
    mol2 = Molecule(smiles="CC(=O)CO")
    mol3 = Molecule(smiles="CCC")
    query = make_stock_query([mol1, mol2, mol3], price={mol1: 10, mol2: 5})
    stock.load(query, "stock1")
    stock.select(["stock1"])

    stock.set_stop_criteria({"price": 5})

    assert mol1 not in stock
    assert mol2 in stock
    assert mol3 in stock


@pytest.fixture
def mocked_mongo_db_query(mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.get_mongo_client")

    def wrapper(**kwargs):
        return mocked_client, MongoDbInchiKeyQuery(**kwargs)

    return wrapper


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


def test_mongodb_integration(stock, mocked_mongo_db_query):
    _, query = mocked_mongo_db_query()
    query.molecules.count_documents.return_value = 1
    stock.load(query, "stock1")
    stock.select(["stock1"])

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
