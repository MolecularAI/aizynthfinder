import logging

from aizynthfinder.aizynthfinder import AiZynthFinder


def state_smiles(state):
    return [mol.smiles for mol in state.mols]


def test_reset_tree():
    finder = AiZynthFinder()
    finder.target_smiles = "CCCO"
    finder.prepare_tree()

    assert finder.tree is not None

    finder.target_smiles = "CCO"

    assert finder.tree is None


def test_one_expansion(setup_aizynthfinder):
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

    nodes = list(finder.tree.graph())
    assert len(nodes) == 2
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert finder.search_stats["iterations"] == 1
    assert finder.search_stats["returned_first"]

    # then test with iteration limit
    finder.config.return_first = False
    finder.config.iteration_limit = 45
    finder.prepare_tree()
    finder.tree_search()

    assert len(finder.tree.graph()) == 2
    assert finder.search_stats["iterations"] == 45
    assert not finder.search_stats["returned_first"]


def test_two_expansions(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
                  |
                child 2
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0},
        child1_smi[1]: {"smiles": ".".join(child2_smi), "prior": 1.0},
    }
    finder = setup_aizynthfinder(lookup, [child1_smi[0], child1_smi[2]] + child2_smi)
    finder.config.return_first = True

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 3
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + child2_smi
    assert finder.search_stats["iterations"] == 1


def test_two_expansions_two_children(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
            /           \
        child 1        child 2
            |             |
        grandchild 1   grandchild 2
    """
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
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2]] + grandchild_smi
    )

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 5
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert (
        state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + grandchild_smi
    )
    assert state_smiles(nodes[3].state) == child2_smi
    assert state_smiles(nodes[4].state) == [child2_smi[0]] + grandchild_smi
    assert finder.search_stats["iterations"] == 100


def test_three_expansions(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
                  |
                child 2
                  |
                child 3 (*)
        - child 3 state is solved
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    child3_smi = ["O=C(Cl)c1ccccc1"]
    lookup = {
        root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0},
        child1_smi[1]: {"smiles": ".".join(child2_smi), "prior": 1.0},
        child2_smi[1]: {"smiles": child3_smi[0], "prior": 1.0},
    }
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2], child2_smi[0]] + child3_smi
    )
    finder.config.return_first = True

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 4
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + child2_smi
    expected_list = [child1_smi[0], child1_smi[2], child2_smi[0]] + child3_smi
    assert state_smiles(nodes[3].state) == expected_list
    assert nodes[3].state.is_solved
    assert finder.search_stats["iterations"] == 1


def test_three_expansions_not_solved(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
                  |
                child 2
                  |
                child 3
        - child 3 state is not solved (not in stock)
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    child3_smi = ["O=C(Cl)c1ccccc1"]
    lookup = {
        root_smi: {"smiles": ".".join(child1_smi), "prior": 1.0},
        child1_smi[1]: {"smiles": ".".join(child2_smi), "prior": 1.0},
        child2_smi[1]: {"smiles": child3_smi[0], "prior": 1.0},
    }
    finder = setup_aizynthfinder(lookup, [child1_smi[0], child1_smi[2], child2_smi[0]])
    finder.config.return_first = True
    finder.config.max_transforms = 2
    finder.config.iteration_limit = 15

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 4
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + child2_smi
    expected_list = [child1_smi[0], child1_smi[2], child2_smi[0]] + child3_smi
    assert state_smiles(nodes[3].state) == expected_list
    assert not nodes[3].state.is_solved
    assert finder.search_stats["iterations"] == 15


def test_two_expansions_no_expandable_root(setup_aizynthfinder):
    """
    Test the following scenario:
                root
                  |
              child 1 (+)

        - child 1 will be selected first for expansion (iteration 1)
        - it has no children that can be expanded (marked by +)
        -- end of iteration 1
        - iteration 2 starts but selecting a leaf will raise an exception
        -- will continue to iterate until reached number of iteration (set 10 in the test)
        * nodes in tree will be root, child 1
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    lookup = {
        root_smi: {
            "smiles": ".".join(child1_smi),
            "prior": 1.0,
        },
        child1_smi[1]: {
            "smiles": "",
            "prior": 0.3,
        },
    }
    finder = setup_aizynthfinder(lookup, [child1_smi[0], child1_smi[2]])
    finder.config.return_first = True
    finder.config.iteration_limit = 10

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 2
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert finder.search_stats["iterations"] == 10


def test_two_expansions_no_reactants_first_child(setup_aizynthfinder):
    """
    Test the following scenario:
                root
            /           \
        child 1 (+)        child 2
                             |
                        grandchild 1 (*)

        - child 1 will be selected first for expansion (iteration 1)
        - it has no children that can be expanded (marked by +)
        -- end of iteration 1
        - child 2 will be selected for expansion  (iteration 2)
        - grandchild 1 will be selected next and it is in stock (marked by *)
        -- a solution is found and the tree search is terminated
        * nodes in tree will be root, child1, child2, grandchild 1
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1CF"]
    grandchild1_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {
            "smiles": "",
            "prior": 0.3,
        },
        child2_smi[1]: {
            "smiles": ".".join(grandchild1_smi),
            "prior": 0.3,
        },
    }
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2]] + grandchild1_smi
    )
    finder.config.return_first = True

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 4
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == child2_smi
    assert state_smiles(nodes[3].state) == [child2_smi[0]] + grandchild1_smi
    assert finder.search_stats["iterations"] == 2


def test_three_expansions_no_reactants_first_child(setup_aizynthfinder):
    """
    Test the following scenario:
                root
            /           \
        child 1 (+)        child 2
                          |
                    grandchild 1
                          |
                    grandchild 2 (*)

        - child 1 will be selected first for expansion (iteration 1)
        - it has no children that can be expanded (marked by +)
        -- end of iteration 1
        - child 2 will be selected for expansion  (iteration 2)
        - grandchild 1 will be selected next
        - grandchild 2 will be selected next and it is in stock (marked by *)
        -- a solution is found and the tree search is terminated
        * nodes in tree will be root, child1, child2, grandchild 1, grandchild 2
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1CF"]
    grandchild1_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    grandchild2_smi = ["O=C(Cl)c1ccccc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {
            "smiles": "",
            "prior": 0.3,
        },
        child2_smi[1]: {
            "smiles": ".".join(grandchild1_smi),
            "prior": 0.3,
        },
        grandchild1_smi[1]: {
            "smiles": ".".join(grandchild2_smi),
            "prior": 1.0,
        },
    }
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2], grandchild1_smi[0]] + grandchild2_smi
    )
    finder.config.return_first = True

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 5
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == child2_smi
    assert state_smiles(nodes[3].state) == [child2_smi[0]] + grandchild1_smi
    expected_list = [child2_smi[0], grandchild1_smi[0]] + grandchild2_smi
    assert state_smiles(nodes[4].state) == expected_list
    assert finder.search_stats["iterations"] == 2


def test_three_expansions_no_reactants_second_level(setup_aizynthfinder):
    """
    Test the following scenario:
                root
            /           \
        child 1         child 2
           |               |
        grandchild 1 (+) grandchild 2 (*)

        - child 1 will be selected first for expansion (iteration 1)
        - grandchild 1 will be selected next,
        - it has no children that can be expanded (marked by x)
        -- end of iteration 1
        - child 2 will be selected for expansion  (iteration 2)
        - grandchild 2 will be selected next and it is in stock (marked by *)
        -- a solution is found and the tree search is terminated
        * nodes in tree will be root, child1, grandchild 1, child2, grandchild 2
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1CF"]
    grandchild1_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    grandchild2_smi = ["N#Cc1cccc(N)c1", "O=C(Cl)c1ccc(F)c(F)c1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {
            "smiles": ".".join(grandchild1_smi),
            "prior": 0.3,
        },
        grandchild1_smi[1]: {
            "smiles": "",
            "prior": 1.0,
        },
        child2_smi[1]: {
            "smiles": ".".join(grandchild2_smi),
            "prior": 0.3,
        },
    }
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2], grandchild1_smi[0]] + grandchild2_smi
    )
    finder.config.return_first = True

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 5
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert (
        state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + grandchild1_smi
    )
    assert state_smiles(nodes[3].state) == child2_smi
    assert state_smiles(nodes[4].state) == [child2_smi[0]] + grandchild2_smi
    assert finder.search_stats["iterations"] == 2


def test_two_expansions_no_reactants_second_child(setup_aizynthfinder):
    """
    Test the following scenario:
                root
            /           \
        child 1        child 2 (+)
            |
        grandchild 1 (*)

        - child 1 will be selected first for expansion (iteration 1)
        - grandchild 1 will be selected next and it is in stock (marked by *)
        -- end of iteration 1
        - child 2 will be selected for expansion  (iteration 2)
        - it has no children that can be expanded (marked with +)
        -- will continue to iterate until reached number of iteration (set 10 in the test)
        * nodes in tree will be root, child1, grandchild 1, child2
    """
    root_smi = "CN1CCC(C(=O)c2cccc(NC(=O)c3ccc(F)cc3)c2F)CC1"
    child1_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F", "O"]
    child2_smi = ["CN1CCC(Cl)CC1", "N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1CF"]
    grandchild1_smi = ["N#Cc1cccc(N)c1F", "O=C(Cl)c1ccc(F)cc1"]
    lookup = {
        root_smi: [
            {"smiles": ".".join(child1_smi), "prior": 0.7},
            {"smiles": ".".join(child2_smi), "prior": 0.3},
        ],
        child1_smi[1]: {
            "smiles": ".".join(grandchild1_smi),
            "prior": 0.3,
        },
        grandchild1_smi[1]: {
            "smiles": "",
            "prior": 1.0,
        },
        child2_smi[1]: {
            "smiles": "",
            "prior": 0.3,
        },
    }
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2]] + grandchild1_smi
    )
    finder.config.iteration_limit = 10

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 4
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert (
        state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + grandchild1_smi
    )
    assert state_smiles(nodes[3].state) == child2_smi
    assert finder.search_stats["iterations"] == 10


def test_two_expansions_cyclic(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
                  |
                child 2
    But making child 2 should be rejected because child 2 == root
    """
    root_smi = "COc1cc2cc(-c3ccc(OC(C)=O)c(OC(C)=O)c3)[n+](C)c(C)c2cc1OC"
    child1_smi = ["COc1cc2cc(-c3ccc(O)c(OC(C)=O)c3)[n+](C)c(C)c2cc1OC"]
    lookup = {
        root_smi: {"smiles": child1_smi[0], "prior": 0.1},
        child1_smi[0]: {
            "smiles": root_smi,
            "prior": 1.0,
        },
    }
    finder = setup_aizynthfinder(lookup, [])
    finder.config.iteration_limit = 1

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 2
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert finder.search_stats["iterations"] == 1


def test_two_expansions_prune_cyclic(setup_aizynthfinder):
    """
    Test the building of this tree:
                root
                  |
                child 1
                  |
                child 2
    Child 2 will not be rejected, but the tree search will not end, so it will
    continue to expand until reaching maximum depth
    """
    root_smi = "COc1cc2cc(-c3ccc(OC(C)=O)c(OC(C)=O)c3)[n+](C)c(C)c2cc1OC"
    child1_smi = ["COc1cc2cc(-c3ccc(O)c(OC(C)=O)c3)[n+](C)c(C)c2cc1OC"]
    lookup = {
        root_smi: {"smiles": child1_smi[0], "prior": 0.1},
        child1_smi[0]: {
            "smiles": root_smi,
            "prior": 1.0,
        },
    }
    finder = setup_aizynthfinder(lookup, [])
    finder.config.iteration_limit = 1
    finder.config.prune_cycles_in_search = False

    finder.tree_search()

    nodes = list(finder.tree.graph())
    assert len(nodes) == 8
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert state_smiles(nodes[2].state) == [root_smi]
    assert finder.search_stats["iterations"] == 1


def test_two_expansions_two_children_one_filtered(setup_aizynthfinder, caplog):
    """
    Test the building of this tree:
                root
            /           \
        child 1        child 2 (*)
            |             |
        grandchild 1   grandchild 2
    child 2 will not be created as that reaction is filtered away
    """
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
    finder = setup_aizynthfinder(
        lookup, [child1_smi[0], child1_smi[2]] + grandchild_smi
    )
    finder.filter_policy[finder.filter_policy.selection[0]].lookup = {
        f"{root_smi}>>{'.'.join(child2_smi)}": 0.2
    }
    finder.config.iteration_limit = 10

    with caplog.at_level(logging.DEBUG):
        finder.tree_search()

    assert not any(
        rec.message.startswith("Reject retro reaction") for rec in caplog.records
    )
    nodes = list(finder.tree.graph())
    assert len(nodes) == 5
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert (
        state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + grandchild_smi
    )
    assert state_smiles(nodes[3].state) == child2_smi
    assert state_smiles(nodes[4].state) == [child2_smi[0]] + grandchild_smi
    assert finder.search_stats["iterations"] == 10

    # Now raise the filter threshold to remove child 2, grandchild 2
    finder.config.filter_cutoff = 0.5
    finder.target_smiles = finder.target_smiles  # Trigger re-set

    with caplog.at_level(logging.DEBUG):
        finder.tree_search()

    assert any(
        rec.message.startswith("Reject retro reaction") for rec in caplog.records
    )
    nodes = list(finder.tree.graph())
    assert len(nodes) == 3
    assert state_smiles(nodes[0].state) == [root_smi]
    assert state_smiles(nodes[1].state) == child1_smi
    assert (
        state_smiles(nodes[2].state) == [child1_smi[0], child1_smi[2]] + grandchild_smi
    )
    assert finder.search_stats["iterations"] == 10
