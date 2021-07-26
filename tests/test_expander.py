from aizynthfinder.aizynthfinder import AiZynthExpander


def test_expander_defaults(get_one_step_expansion, setup_policies):
    expander = AiZynthExpander()
    setup_policies(get_one_step_expansion, config=expander.config)
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"

    reactions = expander.do_expansion(smi)

    assert len(reactions) == 2
    assert len(reactions[0]) == 1
    assert len(reactions[1]) == 1

    assert reactions[0][0].mol.smiles == smi
    assert reactions[1][0].mol.smiles == smi
    assert len(reactions[0][0].reactants[0]) == 2
    assert len(reactions[1][0].reactants[0]) == 2
    smi1 = [mol.smiles for mol in reactions[0][0].reactants[0]]
    smi2 = [mol.smiles for mol in reactions[1][0].reactants[0]]
    assert smi1 != smi2


def test_expander_top1(get_one_step_expansion, setup_policies):
    expander = AiZynthExpander()
    setup_policies(get_one_step_expansion, config=expander.config)
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"

    reactions = expander.do_expansion(smi, return_n=1)

    assert len(reactions) == 1
    smiles_list = [mol.smiles for mol in reactions[0][0].reactants[0]]
    assert smiles_list == ["CCCCOc1ccc(CC(=O)Cl)cc1", "CNO"]


def test_expander_filter(get_one_step_expansion, setup_policies):
    def filter_func(reaction):
        return "CNO" not in [mol.smiles for mol in reaction.reactants[0]]

    expander = AiZynthExpander()
    setup_policies(get_one_step_expansion, config=expander.config)
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"

    reactions = expander.do_expansion(smi, filter_func=filter_func)

    assert len(reactions) == 1
    smiles_list = [mol.smiles for mol in reactions[0][0].reactants[0]]
    assert smiles_list == ["CCCCBr", "CN(O)C(=O)Cc1ccc(O)cc1"]


def test_expander_filter_policy(get_one_step_expansion, setup_policies):
    expander = AiZynthExpander()
    _, filter_strategy = setup_policies(get_one_step_expansion, config=expander.config)
    filter_strategy.lookup[
        "CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"
    ] = 0.5
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"

    reactions = expander.do_expansion(smi)

    assert len(reactions) == 2
    assert reactions[0][0].metadata["feasibility"] == 0.5
    assert reactions[1][0].metadata["feasibility"] == 0.0
