import pytest

from aizynthfinder.search.breadth_first.nodes import MoleculeNode
from aizynthfinder.chem.serialization import MoleculeSerializer, MoleculeDeserializer


@pytest.fixture
def setup_root(default_config):
    def wrapper(smiles):
        return MoleculeNode.create_root(smiles, config=default_config)

    return wrapper


def test_create_root_node(setup_root):
    node = setup_root("CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1")

    assert node.ancestors() == {node.mol}
    assert node.expandable
    assert not node.children


def test_create_stub(setup_root, get_action):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_root(root_smiles)
    reaction = get_action()

    root.add_stub(reaction=reaction)

    assert len(root.children) == 1
    assert len(root.children[0].children) == 2
    rxn_node = root.children[0]
    assert rxn_node.reaction is reaction
    exp_list = [node.mol for node in rxn_node.children]
    assert exp_list == list(reaction.reactants[0])


def test_initialize_stub_one_solved_leaf(
    setup_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, reaction.reactants[0][0])

    root.add_stub(reaction=reaction)

    assert not root.children[0].children[0].expandable
    assert root.children[0].children[1].expandable


def test_serialization_deserialization(
    setup_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, *reaction.reactants[0])
    root.add_stub(reaction=reaction)

    molecule_serializer = MoleculeSerializer()
    dict_ = root.serialize(molecule_serializer)

    molecule_deserializer = MoleculeDeserializer(molecule_serializer.store)
    node = MoleculeNode.from_dict(dict_, default_config, molecule_deserializer)

    assert node.mol == root.mol
    assert len(node.children) == len(root.children)

    rxn_node = node.children[0]
    assert rxn_node.reaction.smarts == reaction.smarts
    assert rxn_node.reaction.metadata == reaction.metadata

    for grandchild1, grandchild2 in zip(rxn_node.children, root.children[0].children):
        assert grandchild1.mol == grandchild2.mol
