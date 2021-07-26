import sys
import gzip
import json
import string
import random

import pytest
import yaml
import pandas as pd

from aizynthfinder.context.config import Configuration
from aizynthfinder.chem import Molecule, TreeMolecule, TemplatedRetroReaction
from aizynthfinder.search.mcts import MctsSearchTree
from aizynthfinder.analysis import TreeAnalysis
from aizynthfinder.search.andor_trees import (
    AndOrSearchTreeBase,
    TreeNodeMixin,
    SplitAndOrTree,
)
from aizynthfinder.context.policy import ExpansionStrategy, FilterStrategy
from aizynthfinder.chem import SmilesBasedRetroReaction
from aizynthfinder.utils.exceptions import RejectionException
from aizynthfinder.chem.serialization import MoleculeDeserializer
from aizynthfinder.aizynthfinder import AiZynthFinder


def pytest_addoption(parser):
    parser.addoption(
        "--run_integration",
        action="store_true",
        default=False,
        help="run integration tests",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run_integration"):
        return

    skip_integration = pytest.mark.skip(reason="need --run_integration option to run")
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(skip_integration)


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "integration: this one is for integration tests."
    )


@pytest.fixture
def add_cli_arguments():
    saved_argv = list(sys.argv)

    def wrapper(args):
        sys.argv = [sys.argv[0]] + args.split(" ")

    yield wrapper
    sys.argv = saved_argv


@pytest.fixture
def create_dummy_smiles_source(tmpdir):
    data = pd.DataFrame(
        {
            "smiles": ["c1ccccc1", "Cc1ccccc1", "c1ccccc1", "CCO"],
            "name": ["benzene", "toluene", "benzene", "ethanol"],
        }
    )

    def wrapper(format_):
        if format_ == "csv":
            filename = str(tmpdir / "smiles_source.csv")
            data.to_csv(filename, index=False)
        else:
            filename = str(tmpdir / "smiles_source.txt")
            with open(filename, "w") as fileobj:
                fileobj.write("\n".join(data["smiles"]))
        return filename

    return wrapper


@pytest.fixture
def create_dummy_stock1(tmpdir):
    data = pd.DataFrame(
        {
            "name": ["benzene", "toluene"],
            "CAS ID": ["71-43-2", "108-88-3"],
            "inchi_key": ["UHOVQNZJYSORNB-UHFFFAOYSA-N", "YXFVVABEGXRONW-UHFFFAOYSA-N"],
        }
    )

    def wrapper(format_):
        if format_ == "hdf5":
            filename = str(tmpdir / "stock1.hdf5")
            data.to_hdf(filename, "table")
        elif format_ == "csv":
            filename = str(tmpdir / "stock1.csv")
            data["inchi_key"].to_csv(filename, index=False)
        else:
            filename = str(tmpdir / "stock1.txt")
            with open(filename, "w") as fileobj:
                fileobj.write("\n".join(data["inchi_key"]))
        return filename

    return wrapper


@pytest.fixture
def create_dummy_stock2(tmpdir):
    filename = str(tmpdir / "stock2.hdf5")
    pd.DataFrame(
        {
            "name": ["benzene", "phenol"],
            "CAS ID": ["71-43-2", "108-95-2"],
            "inchi_key": ["UHOVQNZJYSORNB-UHFFFAOYSA-N", "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"],
        }
    ).to_hdf(filename, "table")
    return filename


@pytest.fixture
def create_dummy_templates(tmpdir):
    def wrapper(ntemplates):
        data = {
            "retro_template": [
                "".join(random.choices(string.ascii_lowercase, k=10))
                for _ in range(ntemplates)
            ]
        }
        filename = str(tmpdir / "dummy_templates.hdf5")
        pd.DataFrame(data).to_hdf(filename, "table")
        return filename

    return wrapper


@pytest.fixture
def default_config():
    return Configuration()


@pytest.fixture
def get_action():
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    mol = TreeMolecule(smiles=smi, parent=None)

    def wrapper(applicable=True):
        if applicable:
            smarts = (
                "([#8:4]-[N;H0;D3;+0:5](-[C;D1;H3:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])"
                ">>(Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([#8:4]-[NH;D2;+0:5]-[C;D1;H3:6])"
            )
        else:
            smarts = (
                "([C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H03])>>"
                "(O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([C:4]-[NH;D2;+0:5]-[C:6])"
            )
        return TemplatedRetroReaction(mol, smarts=smarts, metadata={"dummy": 1})

    return wrapper


@pytest.fixture
def get_branched_expansion():
    return {
        "OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1": {
            "smiles": "OOc1ccc(-c2ccccc2)cc1.NC1CCCC(C2C=CC=C2)C1",
            "prior": 1.0,
        },
        "OOc1ccc(-c2ccccc2)cc1": {"smiles": "c1ccccc1.OOc1ccccc1", "prior": 0.9},
        "OOc1ccccc1": {"smiles": "O.Oc1ccccc1", "prior": 1.0},
        "NC1CCCC(C2C=CC=C2)C1": {"smiles": "NC1CCCCC1.C1=CCC=C1", "prior": 0.1},
    }


@pytest.fixture
def get_expansion_strategy(default_config):
    class LookupExpansionStrategy(ExpansionStrategy):
        _required_kwargs = ["lookup"]

        def __init__(self, key, config, **kwargs):
            super().__init__(key, config, **kwargs)
            self.lookup = kwargs["lookup"]

        def get_actions(self, molecules):
            possible_actions = []
            priors = []
            for mol in molecules:
                if mol.smiles not in self.lookup:
                    continue
                expansion_options = self.lookup[mol.smiles]
                if isinstance(expansion_options, dict):
                    expansion_options = [expansion_options]
                possible_actions.extend(
                    SmilesBasedRetroReaction(mol, reactants_str=item["smiles"])
                    for item in expansion_options
                )
                priors.extend(item["prior"] for item in expansion_options)
            return possible_actions, priors

    def wrapper(lookup, config=None):
        return LookupExpansionStrategy("dummy", config or default_config, lookup=lookup)

    return wrapper


@pytest.fixture
def get_filter_strategy(default_config):
    class LookupFilterStrategy(FilterStrategy):
        _required_kwargs = ["lookup"]

        def __init__(self, key, config, **kwargs):
            super().__init__(key, config, **kwargs)
            self.lookup = kwargs["lookup"]

        def apply(self, reaction):
            prob = self.lookup.get(reaction.smiles)
            if prob is None:
                return
            if prob < self._config.filter_cutoff:
                raise RejectionException(f"Reject {reaction} with prob {prob}")

        def feasibility(self, reaction):
            prob = self.lookup.get(reaction.smiles)
            if prob is None:
                return False, 0.0
            return prob < self._config.filter_cutoff, prob

    def wrapper(lookup, config=None):
        return LookupFilterStrategy("dummy", config or default_config, lookup=lookup)

    return wrapper


@pytest.fixture
def get_linear_expansion():
    return {
        "OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1": {
            "smiles": "OOc1ccc(-c2ccccc2)cc1.NC1CCCC(C2C=CC=C2)C1",
            "prior": 1.0,
        },
        "OOc1ccc(-c2ccccc2)cc1": {"smiles": "c1ccccc1.OOc1ccccc1", "prior": 1.0},
    }


@pytest.fixture
def get_one_step_expansion():
    return {
        "CCCCOc1ccc(CC(=O)N(C)O)cc1": [
            {"smiles": "CCCCOc1ccc(CC(=O)Cl)cc1.CNO", "prior": 0.7},
            {"smiles": "CCCCBr.CN(O)C(=O)Cc1ccc(O)cc1", "prior": 0.5},
            {"smiles": "", "prior": 0.3},
        ]
    }


@pytest.fixture
def load_reaction_tree(shared_datadir):
    def wrapper(filename, index=0):
        filename = str(shared_datadir / filename)
        with open(filename, "r") as fileobj:
            trees = json.load(fileobj)
        if isinstance(trees, dict):
            return trees
        else:
            return trees[index]

    return wrapper


@pytest.fixture
def load_mcts_tree(default_config, shared_datadir):
    def wrapper(filename, config=None):
        filename = str(shared_datadir / filename)
        return MctsSearchTree.from_json(filename, config or default_config)

    return wrapper


@pytest.fixture
def setup_aizynthfinder(setup_policies, setup_stock):
    def wrapper(expansions, stock):
        finder = AiZynthFinder()
        root_smi = list(expansions.keys())[0]
        setup_policies(expansions, config=finder.config)
        setup_stock(finder.config, *stock)
        finder.target_smiles = root_smi
        return finder

    return wrapper


@pytest.fixture
def setup_analysis(default_config, shared_datadir, tmpdir, setup_stock):
    setup_stock(
        default_config, "N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1", "CN1CCC(Cl)CC1", "O"
    )

    with gzip.open(shared_datadir / "full_search_tree.json.gz", "rb") as gzip_obj:
        with open(tmpdir / "full_search_tree.json", "wb") as fileobj:
            fileobj.write(gzip_obj.read())
    tree = MctsSearchTree.from_json(tmpdir / "full_search_tree.json", default_config)
    nodes = list(tree.graph())

    def wrapper(scorer=None):
        return TreeAnalysis(tree, scorer=scorer), nodes

    return wrapper


@pytest.fixture
def setup_analysis_andor_tree(default_config, shared_datadir, setup_stock):  # noqa
    setup_stock(
        default_config,
        "Nc1ccc(NC(=S)Nc2ccccc2)cc1",
        "Cc1ccc2nc3ccccc3c(Cl)c2c1",
        "Nc1ccccc1",
        "Nc1ccc(N=C=S)cc1",
        "Cc1ccc2nc3ccccc3c(Br)c2c1",
        "Nc1ccc(Br)cc1",
    )

    class BasicAndOrTree(AndOrSearchTreeBase):
        def __init__(self, filename, config):
            super().__init__(config)
            self._mol_nodes = []
            with open(filename, "r") as fileobj:
                dict_ = json.load(fileobj)
            mol_deser = MoleculeDeserializer(dict_["molecules"])
            self.root = AndOrNode(dict_["tree"], config, mol_deser, self)

        @property
        def mol_nodes(self):
            return self._mol_nodes

        def one_iteration(self):
            return False

        def routes(self):
            return SplitAndOrTree(self.root, self.config.stock).routes

    class AndOrNode(TreeNodeMixin):
        def __init__(self, dict_, config, molecules, tree):
            self.tree = tree
            self.config = config
            self._children = [
                AndOrNode(child, config, molecules, tree) for child in dict_["children"]
            ]
            if "reaction" in dict_:
                self._obj = TemplatedRetroReaction(
                    molecules[dict_["reaction"]["mol"]],
                    smarts=dict_["reaction"]["smarts"],
                    index=dict_["reaction"]["index"],
                    metadata=dict_["reaction"].get("metadata", {}),
                )
                self._solved = all(child.prop["solved"] for child in self._children)
            else:
                self.tree._mol_nodes.append(self)
                self._obj = molecules[dict_["mol"]]
                self._solved = self._obj in self.config.stock

        @property
        def prop(self):
            obj_key = (
                "reaction" if isinstance(self._obj, TemplatedRetroReaction) else "mol"
            )
            return {obj_key: self._obj, "solved": self._solved}

        @property
        def children(self):
            return self._children

    tree = BasicAndOrTree(str(shared_datadir / "and_or_tree.json"), default_config)

    def wrapper(scorer=None):
        return TreeAnalysis(tree, scorer=scorer)

    return wrapper


@pytest.fixture
def setup_expanded_mcts(default_config, setup_policies):
    def wrapper(expansions):
        expansion_strategy, _ = setup_policies(expansions)
        root_smiles = list(expansion_strategy.lookup.keys())[0]
        tree = MctsSearchTree(config=default_config, root_smiles=root_smiles)
        tree.one_iteration()
        nodes = [node for node in tree.graph()]
        return tree, nodes[-1]

    return wrapper


@pytest.fixture
def setup_branched_mcts(setup_stock, setup_expanded_mcts, get_branched_expansion):
    def wrapper(exclude_from_stock=None):
        exclude_from_stock = exclude_from_stock or []
        stock = [
            smi
            for smi in ["c1ccccc1", "O", "Oc1ccccc1", "NC1CCCCC1", "C1=CCC=C1"]
            if smi not in exclude_from_stock
        ]
        setup_stock(None, *stock)
        for smi in exclude_from_stock:
            get_branched_expansion[smi] = {"smiles": "", "prior": 1}
        return setup_expanded_mcts(get_branched_expansion)

    return wrapper


@pytest.fixture
def setup_branched_reaction_tree(setup_branched_mcts):
    def wrapper(exclude_from_stock=None):
        _, node = setup_branched_mcts(exclude_from_stock or [])
        return node.to_reaction_tree()

    return wrapper


@pytest.fixture
def setup_linear_mcts(setup_stock, setup_expanded_mcts, get_linear_expansion):
    def wrapper(exclude_from_stock=None):
        exclude_from_stock = exclude_from_stock or []
        stock = [
            smi
            for smi in ["NC1CCCC(C2C=CC=C2)C1", "c1ccccc1", "OOc1ccccc1"]
            if smi not in exclude_from_stock
        ]
        setup_stock(None, *stock)
        for smi in exclude_from_stock:
            get_linear_expansion[smi] = {"smiles": "", "prior": 1}
        return setup_expanded_mcts(get_linear_expansion)

    return wrapper


@pytest.fixture
def setup_linear_reaction_tree(setup_linear_mcts):
    def wrapper(exclude_from_stock=None):
        _, node = setup_linear_mcts(exclude_from_stock or [])
        return node.to_reaction_tree()

    return wrapper


@pytest.fixture
def setup_policies(default_config, get_filter_strategy, get_expansion_strategy):
    def wrapper(expansions, filters=None, config=None):
        config = config or default_config

        expansion_strategy = get_expansion_strategy(expansions, config)
        config.expansion_policy.load(expansion_strategy)
        config.expansion_policy.select("dummy")

        filter_strategy = get_filter_strategy(filters or {}, config=config)
        config.filter_policy.load(filter_strategy)
        config.filter_policy.select("dummy")
        return expansion_strategy, filter_strategy

    return wrapper


@pytest.fixture
def setup_stock(default_config, tmpdir):
    """
    Fixture for setting up stock of inchi keys in a textfile.
    Will return a function that should be called with any number of Molecule objects as arguments
    """

    def wrapper(config=None, *molecules):
        config = config or default_config
        molecules = [
            Molecule(smiles=mol) if isinstance(mol, str) else mol for mol in molecules
        ]
        filename = str(tmpdir / "stock.txt")
        with open(filename, "w") as fileobj:
            fileobj.write("\n".join([mol.inchi_key for mol in molecules]))
        config.stock.load(filename, "stock")
        config.stock.select("stock")

    return wrapper


@pytest.fixture
def write_yaml(tmpdir):
    filename = str(tmpdir / "test.yaml")

    def wrapper(dict_):
        with open(filename, "w") as fileobj:
            yaml.dump(dict_, fileobj)
        return filename

    return wrapper
