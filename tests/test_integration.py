import json
import random
import os

import pytest

from aizynthfinder.aizynthfinder import AiZynthFinder


def _remove_meta_data(dict_):
    if "metadata" in dict_:
        del dict_["metadata"]
    if "children" in dict_:
        for child in dict_["children"]:
            _remove_meta_data(child)


@pytest.fixture(scope="module")
def read_options(request):
    policy_model = os.environ.get("TEST_POLICY_MODEL")
    policy_templates = os.environ.get("TEST_POLICY_TEMPLATES")
    stock_file = os.environ.get("TEST_STOCK")

    if not (policy_model and policy_templates and stock_file):
        pytest.xfail(
            "You need to specify the environmental variables TEST_POLICY_MODEL\
                and TEST_POLICY_TEMPLATES and TEST_STOCK to run integration tests"
        )

    return policy_model, policy_templates, stock_file


@pytest.fixture(scope="module")
def setup_finder(read_options):
    policy_model, policy_templates, stock_file = read_options
    random.seed(350)
    finder = AiZynthFinder()
    finder.config.stock.load_stock(stock_file, "test_stock")
    finder.config.stock.select_stocks("test_stock")
    finder.config.policy.load_policy(policy_model, policy_templates, "test")
    finder.config.policy.select_policy("test")
    return finder


@pytest.fixture
def finder_output(shared_datadir):
    filename = str(shared_datadir / "finder_output.json")
    with open(filename, "r") as fileobj:
        expected_outputs = json.load(fileobj)

    def wrapper(smiles):
        expected_output = expected_outputs[smiles]

        trees_filename = str(shared_datadir / expected_output["trees"])
        with open(trees_filename, "r") as fileobj:
            expected_output["trees"] = json.load(fileobj)

        return expected_output

    return wrapper


@pytest.mark.integration
@pytest.mark.parametrize(
    "smiles",
    [
        "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1",
        "COc1cc2cc(-c3ccc(O)c(O)c3)[n+](C)c(C)c2cc1OC",
    ],
)
def test_full_flow(finder_output, setup_finder, smiles):
    expected_output = finder_output(smiles)
    finder = setup_finder
    finder.target_smiles = smiles
    finder.prepare_tree()
    finder.tree_search()
    finder.build_routes()

    best_score = finder.routes.scores
    assert all(
        abs(actual - expected) < 0.001
        for actual, expected in zip(best_score, expected_output["best_scores"])
    )

    trees = finder.routes.dicts
    for tree in trees:
        _remove_meta_data(tree)
    assert all(tree in expected_output["trees"] for tree in trees)
    assert all(tree in trees for tree in expected_output["trees"])


@pytest.mark.integration
def test_run_from_json(finder_output, setup_finder):
    smiles = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    params = {
        "stocks": ["test_stock"],
        "policy": "test",
        "C": setup_finder.config.C,
        "max_transforms": setup_finder.config.max_transforms,
        "cutoff_cumulative": setup_finder.config.cutoff_cumulative,
        "cutoff_number": setup_finder.config.cutoff_number,
        "smiles": smiles,
        "return_first": setup_finder.config.return_first,
        "time_limit": setup_finder.config.time_limit,
        "iteration_limit": setup_finder.config.iteration_limit,
        "exclude_target_from_stock": True,
    }
    expected_output = finder_output(smiles)

    result = setup_finder.run_from_json(params)

    for tree in result["trees"]:
        _remove_meta_data(tree)
    assert result["request"] == params
    assert all(tree in expected_output["trees"] for tree in result["trees"])
    assert all(tree in result["trees"] for tree in expected_output["trees"])
