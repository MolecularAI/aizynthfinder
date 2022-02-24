import random

import pytest

from aizynthfinder.search.breadth_first.search_tree import SearchTree


def test_one_iteration(default_config, setup_policies, setup_stock):
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

    assert len(tree.mol_nodes) == 1

    assert not tree.one_iteration()

    assert len(tree.mol_nodes) == 6
    smiles = [node.mol.smiles for node in tree.mol_nodes]
    assert smiles == [root_smi] + child1_smi + child2_smi

    assert tree.one_iteration()

    assert len(tree.mol_nodes) == 10
    smiles = [node.mol.smiles for node in tree.mol_nodes]
    assert (
        smiles == [root_smi] + child1_smi + child2_smi + grandchild_smi + grandchild_smi
    )

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

    assert len(tree.mol_nodes) == 1

    tree.one_iteration()
    assert len(tree.mol_nodes) == 6

    assert not tree.one_iteration()

    assert len(tree.mol_nodes) == 8

    with pytest.raises(StopIteration):
        tree.one_iteration()


def test_routes(default_config, setup_policies, setup_stock):
    random.seed(666)
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["O", "CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F"]
    grandchild_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {"smiles": ".".join(grandchild_smi), "prior": 0.7},
    }
    stock = [child1_smi[0], child1_smi[2]] + grandchild_smi
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
    smiles = [mol.smiles for mol in routes[1].molecules()]
    assert smiles == [root_smi] + child1_smi + grandchild_smi
    smiles = [mol.smiles for mol in routes[0].molecules()]
    assert smiles == [root_smi] + child2_smi + grandchild_smi
