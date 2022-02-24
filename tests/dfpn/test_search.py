import random

import pytest

from aizynthfinder.search.dfpn.search_tree import SearchTree


def test_search(default_config, setup_policies, setup_stock):
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F"]
    grandchild_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {"smiles": ".".join(grandchild_smi), "prior": 0.7},
        child2_smi[1]: {"smiles": ".".join(grandchild_smi), "prior": 0.7},
    }
    stock = [child1_smi[0], child1_smi[2]] + grandchild_smi
    setup_policies(lookup, config=default_config)
    setup_stock(default_config, *stock)
    tree = SearchTree(default_config, root_smi)

    assert not tree.one_iteration()

    routes = tree.routes()
    assert all([not route.is_solved for route in routes])
    assert len(tree.mol_nodes) == 4

    assert not tree.one_iteration()
    assert tree.one_iteration()
    assert tree.one_iteration()
    assert tree.one_iteration()
    assert tree.one_iteration()

    routes = tree.routes()
    assert len(routes) == 2
    assert all(route.is_solved for route in routes)

    with pytest.raises(StopIteration):
        tree.one_iteration()


def test_search_incomplete(default_config, setup_policies, setup_stock):
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["ClC(=O)c1ccc(F)cc1", "CN1CCC(CC1)C(=O)c1cccc(N)c1F"]
    grandchild_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {"smiles": ".".join(grandchild_smi), "prior": 0.7},
    }
    stock = [child1_smi[0], child1_smi[2]] + [grandchild_smi[0]]
    setup_policies(lookup, config=default_config)
    setup_stock(default_config, *stock)
    tree = SearchTree(default_config, root_smi)

    while True:
        try:
            tree.one_iteration()
        except StopIteration:
            break

    routes = tree.routes()
    assert len(routes) == 2
    assert all(not route.is_solved for route in routes)
