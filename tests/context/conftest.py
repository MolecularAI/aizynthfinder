from typing import Dict, List

import numpy as np
import pandas as pd
import pytest
import pytest_mock

from aizynthfinder.context.policy import TemplateBasedExpansionStrategy
from aizynthfinder.context.scoring import BrokenBondsScorer
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
    def wrapper(
        key="policy1", templates=None, expansion_cls=TemplateBasedExpansionStrategy
    ):
        if templates is None:
            templates_filename = create_dummy_templates(3)
        else:
            templates_filename = create_templates_file(templates)

        strategy = expansion_cls(
            key, default_config, model="dummy.onnx", template=templates_filename
        )

        return strategy, mock_onnx_model

    return wrapper


@pytest.fixture
def setup_mcts_broken_bonds(setup_stock, setup_expanded_mcts, shared_datadir):
    def wrapper(broken=True, config=None):
        root_smi = "CN1CCC(C(=O)c2cccc([NH:1][C:2](=O)c3ccc(F)cc3)c2F)CC1"

        reaction_template = pd.read_csv(
            shared_datadir / "test_reactions_template.csv", sep="\t"
        )
        template1_smarts = reaction_template["RetroTemplate"][0]
        template2_smarts = reaction_template["RetroTemplate"][1]

        child1_smi = ["N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "CN1CCC(Cl)CC1", "O"]
        child2_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]

        if not broken:
            template2_smarts = reaction_template["RetroTemplate"][2]
            child2_smi = ["N#Cc1cccc(Cl)c1F", "NC(=O)c1ccc(F)cc1"]

        lookup = {
            root_smi: {"smarts": template1_smarts, "prior": 1.0},
            child1_smi[0]: {
                "smarts": template2_smarts,
                "prior": 1.0,
            },
        }

        stock = [child1_smi[1], child1_smi[2]] + child2_smi

        if config:
            config.scorers.create_default_scorers()
            config.scorers.load(BrokenBondsScorer(config))
        setup_stock(config, *stock)
        return setup_expanded_mcts(lookup)

    return wrapper
