import numpy as np
import networkx as nx

from aizynthfinder.search.retrostar.nodes import MoleculeNode
from aizynthfinder.chem.serialization import MoleculeSerializer, MoleculeDeserializer
from aizynthfinder.search.andor_trees import ReactionTreeFromAndOrTrace


def test_create_root_node(setup_star_root):
    node = setup_star_root("CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1")

    assert node.target_value == 0
    assert node.ancestors() == {node.mol}
    assert node.expandable
    assert not node.children


def test_close_single_node(setup_star_root):
    node = setup_star_root("CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1")

    assert node.expandable

    delta = node.close()

    assert not np.isfinite(delta)
    assert not node.solved
    assert not node.expandable


def test_create_stub(setup_star_root, get_action):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_star_root(root_smiles)
    reaction = get_action()

    root.add_stub(cost=5.0, reaction=reaction)

    assert len(root.children) == 1
    assert len(root.children[0].children) == 2
    rxn_node = root.children[0]
    assert rxn_node.reaction is reaction
    exp_list = [node.mol for node in rxn_node.children]
    assert exp_list == list(reaction.reactants[0])
    assert rxn_node.value == rxn_node.target_value == 5
    assert not rxn_node.solved

    # This is done after a node has been expanded
    delta = root.close()

    assert delta == 5.0
    assert root.value == 5.0
    assert rxn_node.children[0].ancestors() == {root.mol, rxn_node.children[0].mol}


def test_initialize_stub_one_solved_leaf(
    setup_star_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_star_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, reaction.reactants[0][0])

    root.add_stub(cost=5.0, reaction=reaction)
    root.close()

    assert not root.children[0].solved
    assert not root.solved
    assert root.children[0].children[0].solved


def test_initialize_stub_two_solved_leafs(
    setup_star_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_star_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, *reaction.reactants[0])

    root.add_stub(cost=5.0, reaction=reaction)
    root.close()

    assert root.children[0].solved
    assert root.solved


def test_serialization_deserialization(
    setup_star_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_star_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, *reaction.reactants[0])
    root.add_stub(cost=5.0, reaction=reaction)
    root.close()

    molecule_serializer = MoleculeSerializer()
    dict_ = root.serialize(molecule_serializer)

    molecule_deserializer = MoleculeDeserializer(molecule_serializer.store)
    node = MoleculeNode.from_dict(dict_, default_config, molecule_deserializer)

    assert node.mol == root.mol
    assert node.value == root.value
    assert node.cost == root.cost
    assert len(node.children) == len(root.children)

    rxn_node = node.children[0]
    assert rxn_node.reaction.smarts == reaction.smarts
    assert rxn_node.reaction.metadata == reaction.metadata
    assert rxn_node.cost == root.children[0].cost
    assert rxn_node.value == root.children[0].value

    for grandchild1, grandchild2 in zip(rxn_node.children, root.children[0].children):
        assert grandchild1.mol == grandchild2.mol


def test_converstion_to_reaction_tree(
    setup_star_root, get_action, default_config, setup_stock
):
    root_smiles = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    root = setup_star_root(root_smiles)
    reaction = get_action()
    setup_stock(default_config, *reaction.reactants[0])
    root.add_stub(cost=5.0, reaction=reaction)
    root.close()
    graph = nx.DiGraph()
    graph.add_edge(root, root.children[0])
    graph.add_edge(root.children[0], root.children[0].children[0])
    graph.add_edge(root.children[0], root.children[0].children[1])

    rt = ReactionTreeFromAndOrTrace(graph, default_config.stock).tree

    molecules = list(rt.molecules())
    rt_reactions = list(rt.reactions())
    assert len(molecules) == 3
    assert len(list(rt.leafs())) == 2
    assert len(rt_reactions) == 1
    assert molecules[0].inchi_key == root.mol.inchi_key
    assert molecules[1].inchi_key == root.children[0].children[0].mol.inchi_key
    assert molecules[2].inchi_key == root.children[0].children[1].mol.inchi_key
    assert rt_reactions[0].reaction_smiles() == reaction.reaction_smiles()
