from aizynthfinder.search.retrostar.search_tree import SearchTree
from aizynthfinder.chem.serialization import MoleculeSerializer


def test_one_iteration(setup_search_tree):
    tree = setup_search_tree

    tree.one_iteration()

    assert len(tree.root.children) == 1
    assert len(tree.root.children[0].children) == 3


def test_one_iteration_filter_unfeasible(setup_search_tree):
    tree = setup_search_tree
    smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1>>CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O"
    tree.config.filter_policy["dummy"].lookup[smi] = 0.0

    tree.one_iteration()
    assert len(tree.root.children) == 0


def test_one_iteration_filter_feasible(setup_search_tree):
    tree = setup_search_tree
    smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1>>CN1CCC(Cl)CC1.N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F.O"
    tree.config.filter_policy["dummy"].lookup[smi] = 0.5

    tree.one_iteration()
    assert len(tree.root.children) == 1


def test_one_expansion_with_finder(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    lookup = {root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0}}
    finder = setup_aizynthfinder(lookup, child1_smi)

    # Test first with return_first
    finder.config.return_first = True
    finder.tree_search()

    nodes = finder.tree.mol_nodes
    assert len(nodes) == 4
    assert nodes[0].mol.smiles == root_smi
    assert nodes[1].mol.smiles == child1_smi[0]
    assert finder.search_stats["iterations"] == 1
    assert finder.search_stats["returned_first"]

    # then test with iteration limit
    finder.config.return_first = False
    finder.config.iteration_limit = 45
    finder.prepare_tree()
    finder.tree_search()

    assert len(finder.tree.mol_nodes) == 4
    # It will not continue because it cannot expand any more nodes
    assert finder.search_stats["iterations"] == 2
    assert not finder.search_stats["returned_first"]


def test_serialization_deserialization(
    mocker, setup_search_tree, tmpdir, default_config
):
    tree = setup_search_tree
    tree.one_iteration()

    mocked_json_dump = mocker.patch(
        "aizynthfinder.search.retrostar.search_tree.json.dump"
    )
    serializer = MoleculeSerializer()
    filename = str(tmpdir / "dummy.json")

    # Test serialization

    tree.serialize(filename)

    expected_dict = {
        "tree": tree.root.serialize(serializer),
        "molecules": serializer.store,
    }

    mocked_json_dump.assert_called_once_with(
        expected_dict, mocker.ANY, indent=mocker.ANY
    )

    # Test deserialization

    mocker.patch(
        "aizynthfinder.search.retrostar.search_tree.json.load",
        return_value=expected_dict,
    )
    mocker.patch(
        "aizynthfinder.search.retrostar.nodes.deserialize_action", return_value=None
    )

    new_tree = SearchTree.from_json(filename, default_config)

    assert new_tree.root.mol == tree.root.mol
    assert len(new_tree.root.children) == len(tree.root.children)


def test_split_andor_tree(shared_datadir, default_config):
    tree = SearchTree.from_json(
        str(shared_datadir / "andor_tree_for_clustering.json"), default_config
    )

    routes = tree.routes()

    assert len(routes) == 3


def test_update(shared_datadir, default_config, setup_stock):
    setup_stock(
        default_config,
        "Nc1ccc(NC(=S)Nc2ccccc2)cc1",
        "Cc1ccc2nc3ccccc3c(Cl)c2c1",
        "Nc1ccccc1",
        "Nc1ccc(N=C=S)cc1",
        "Cc1ccc2nc3ccccc3c(Br)c2c1",
        "Nc1ccc(Br)cc1",
    )
    tree = SearchTree.from_json(
        str(shared_datadir / "andor_tree_for_clustering.json"), default_config
    )

    saved_root_value = tree.root.value
    tree.mol_nodes[-1].parent.update(35, from_mol=tree.mol_nodes[-1].mol)

    assert [child.value for child in tree.root.children] == [5.0, 45.0]
    assert tree.root.value != saved_root_value
    assert tree.root.value == 5

    tree.serialize("temp.json")
