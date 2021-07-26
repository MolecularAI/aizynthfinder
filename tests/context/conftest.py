import pytest
import numpy as np

from aizynthfinder.context.policy import TemplateBasedExpansionStrategy
from aizynthfinder.context.stock import (
    MongoDbInchiKeyQuery,
    StockQueryMixin,
    StockException,
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


@pytest.fixture
def mock_keras_model(mocker):
    class MockedKerasModel(mocker.MagicMock):
        @property
        def input(self):
            pass

        @property
        def output(self):
            pass

        def predict(self, *_):
            pass

    mocker.patch.object(
        MockedKerasModel, "input", mocker.PropertyMock(return_value=np.zeros((3, 3)))
    )
    mocker.patch.object(
        MockedKerasModel, "output", mocker.PropertyMock(return_value=np.zeros((3, 3)))
    )
    mocker.patch.object(
        MockedKerasModel,
        "predict",
        mocker.MagicMock(return_value=np.array([[0.2, 0.7, 0.1]])),
    )

    return mocker.patch(
        "aizynthfinder.utils.models.load_keras_model", return_value=MockedKerasModel
    )


@pytest.fixture
def mocked_mongo_db_query(mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.queries.get_mongo_client")

    def wrapper(**kwargs):
        return mocked_client, MongoDbInchiKeyQuery(**kwargs)

    return wrapper


@pytest.fixture
def setup_template_expansion_policy(
    default_config, create_dummy_templates, mock_keras_model
):
    templates_filename = create_dummy_templates(3)

    def wrapper(key="policy1"):
        strategy = TemplateBasedExpansionStrategy(
            key, default_config, source="dummy.hdf5", templatefile=templates_filename
        )

        return strategy, mock_keras_model

    return wrapper


@pytest.fixture
def setup_stock_with_query(default_config, create_dummy_stock1):
    stock = default_config.stock

    def wrapper(query=None):
        query = query if query is not None else create_dummy_stock1("hdf5")
        stock.load(query, "stock1")
        stock.select(["stock1"])
        return stock

    return wrapper
