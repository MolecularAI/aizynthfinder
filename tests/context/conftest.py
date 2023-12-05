from typing import Dict, List

import numpy as np
import pandas as pd
import pytest
import pytest_mock

from aizynthfinder.context.policy import TemplateBasedExpansionStrategy
from aizynthfinder.context.stock import (
    MongoDbInchiKeyQuery,
    StockException,
    StockQueryMixin,
)


@pytest.fixture
def create_templates_file(tmpdir):
    def wrapper(templates):
        data = {"retro_template": templates}
        filename = str(tmpdir / "dummy_templates.hdf5")
        pd.DataFrame(data).to_hdf(filename, "table")
        return filename

    return wrapper


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
def mocked_mongo_db_query(mocker):
    mocked_client = mocker.patch("aizynthfinder.context.stock.queries.get_mongo_client")

    def wrapper(**kwargs):
        return mocked_client, MongoDbInchiKeyQuery(**kwargs)

    return wrapper


@pytest.fixture
def setup_template_expansion_policy(
    default_config, create_dummy_templates, create_templates_file, mock_onnx_model
):
    def wrapper(key="policy1", templates=None):
        if templates is None:
            templates_filename = create_dummy_templates(3)
        else:
            templates_filename = create_templates_file(templates)

        strategy = TemplateBasedExpansionStrategy(
            key, default_config, model="dummy.onnx", template=templates_filename
        )

        return strategy, mock_onnx_model

    return wrapper
