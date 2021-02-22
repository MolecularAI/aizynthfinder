import json
import random
import os

import pytest

from aizynthfinder.aizynthfinder import AiZynthFinder


def _remove_and_yield_meta_data(dict_):
    if "metadata" in dict_:
        yield dict_["metadata"]
        del dict_["metadata"]
    if "children" in dict_:
        for child in dict_["children"]:
            yield from _remove_and_yield_meta_data(child)


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
    finder.config.stock.load(stock_file, "test_stock")
    finder.config.stock.select("test_stock")
    finder.config.expansion_policy.load(policy_model, policy_templates, "dummy")
    finder.config.expansion_policy.select("dummy")
    return finder


@pytest.fixture
def finder_output(shared_datadir):
    filename = str(shared_datadir / "finder_output.json")
    with open(filename, "r") as fileobj:
        expected_outputs = json.load(fileobj)

    def wrapper(smiles):
        expected_output = expected_outputs[smiles]
        return expected_output

    return wrapper


@pytest.mark.integration
@pytest.mark.parametrize(
    "smiles",
    [
        "CCn1nc(CC(C)C)cc1C(=O)NCc1c(C)cc(C)nc1OC",
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

    assert len(finder.routes) == len(expected_output["trees"])

    for actual, expected in zip(finder.routes.scores, expected_output["scores"]):
        assert actual == expected

    for actual, expected in zip(finder.routes.dicts, expected_output["trees"]):
        metalist_actual = [metadata for metadata in _remove_and_yield_meta_data(actual)]
        metalist_expected = [
            metadata for metadata in _remove_and_yield_meta_data(expected)
        ]
        for meta_actual, meta_expected in zip(metalist_actual, metalist_expected):
            assert meta_actual == pytest.approx(meta_expected)
        assert actual == expected


@pytest.mark.integration
@pytest.mark.parametrize(
    "include_scores",
    [False, True],
)
def test_run_from_json(finder_output, setup_finder, include_scores):

    smiles = "CCn1nc(CC(C)C)cc1C(=O)NCc1c(C)cc(C)nc1OC"
    params = {
        "stocks": ["test_stock"],
        "policy": "dummy",
        "C": setup_finder.config.C,
        "max_transforms": setup_finder.config.max_transforms,
        "cutoff_cumulative": setup_finder.config.cutoff_cumulative,
        "cutoff_number": setup_finder.config.cutoff_number,
        "smiles": smiles,
        "return_first": setup_finder.config.return_first,
        "time_limit": setup_finder.config.time_limit,
        "iteration_limit": setup_finder.config.iteration_limit,
        "exclude_target_from_stock": True,
        "filter_cutoff": setup_finder.config.filter_cutoff,
    }
    if include_scores:
        params["score_trees"] = True
    expected_output = finder_output(smiles)

    result = setup_finder.run_from_json(params)

    if include_scores:
        del params["score_trees"]
        for tree, expected in zip(result["trees"], expected_output["all_scores"]):
            assert tree["scores"] == expected
            del tree["scores"]

    assert result["request"] == params

    for actual, expected in zip(result["trees"], expected_output["trees"]):
        metalist_actual = [metadata for metadata in _remove_and_yield_meta_data(actual)]
        metalist_expected = [
            metadata for metadata in _remove_and_yield_meta_data(expected)
        ]
        for meta_actual, meta_expected in zip(metalist_actual, metalist_expected):
            assert meta_actual == pytest.approx(meta_expected)
        assert actual == expected
