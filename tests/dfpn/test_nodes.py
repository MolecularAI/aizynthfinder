import pytest

from aizynthfinder.search.dfpn.nodes import MoleculeNode, BIG_INT
from aizynthfinder.search.dfpn import SearchTree


@pytest.fixture
def setup_root(default_config):
    def wrapper(smiles):
        owner = SearchTree(default_config)
        return MoleculeNode.create_root(smiles, config=default_config, owner=owner)

    return wrapper


def test_create_root_node(setup_root):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")

    assert node.expandable
    assert not node.children
    assert node.dn == 1
    assert node.pn == 1


def test_expand_mol_node(
    default_config, setup_root, setup_policies, get_linear_expansion
):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")
    setup_policies(get_linear_expansion)

    node.expand()

    assert not node.expandable
    assert len(node.children) == 1


def test_promising_child(
    default_config, setup_root, setup_policies, get_linear_expansion
):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")
    setup_policies(get_linear_expansion)
    node.expand()

    child = node.promising_child()

    assert child is node.children[0]
    assert child.pn_threshold == BIG_INT - 1
    assert child.dn_threshold == BIG_INT


def test_expand_reaction_node(
    default_config, setup_root, setup_policies, get_linear_expansion
):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")
    setup_policies(get_linear_expansion)
    node.expand()
    child = node.promising_child()

    child.expand()

    assert len(child.children) == 2


def test_promising_child_reaction_node(
    default_config,
    setup_root,
    setup_policies,
    get_linear_expansion,
):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")
    setup_policies(get_linear_expansion)
    node.expand()

    child = node.promising_child()
    child.expand()

    grandchild = child.promising_child()

    assert grandchild.mol.smiles == "OOc1ccc(-c2ccccc2)cc1"
    assert grandchild.pn_threshold == BIG_INT - 1
    assert grandchild.dn_threshold == 2


def test_update(
    default_config,
    setup_root,
    setup_policies,
    get_linear_expansion,
    setup_stock,
):
    node = setup_root("OOc1ccc(-c2ccc(NC3CCCC(C4C=CC=C4)C3)cc2)cc1")
    setup_stock(default_config, "OOc1ccc(-c2ccccc2)cc1", "NC1CCCC(C2C=CC=C2)C1")
    setup_policies(get_linear_expansion)
    node.expand()

    child = node.promising_child()
    child.expand()
    child.update()

    node.update()

    assert node.proven
    assert not node.disproven
