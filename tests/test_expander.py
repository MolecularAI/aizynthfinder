from aizynthfinder.chem import TreeMolecule
from aizynthfinder.aizynthfinder import AiZynthExpander


def test_expander_defaults(mock_expansion_policy):
    expander = AiZynthExpander()
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    mock_expansion_policy(TreeMolecule(parent=None, smiles=smi))

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


def test_expander_top1(mock_expansion_policy):
    expander = AiZynthExpander()
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    mock_expansion_policy(TreeMolecule(parent=None, smiles=smi))

    reactions = expander.do_expansion(smi, return_n=1)

    assert len(reactions) == 1
    smiles_list = [mol.smiles for mol in reactions[0][0].reactants[0]]
    assert smiles_list == ["CCCCOc1ccc(CC(=O)Cl)cc1", "CNO"]


def test_expander_filter(mock_expansion_policy):
    def filter_func(reaction):
        return "CNO" not in [mol.smiles for mol in reaction.reactants[0]]

    expander = AiZynthExpander()
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    mock_expansion_policy(TreeMolecule(parent=None, smiles=smi))

    reactions = expander.do_expansion(smi, filter_func=filter_func)

    assert len(reactions) == 1
    smiles_list = [mol.smiles for mol in reactions[0][0].reactants[0]]
    assert smiles_list == ["CCCCBr", "CN(O)C(=O)Cc1ccc(O)cc1"]


def test_expander_filter_policy(mock_expansion_policy, mock_policy_model):
    expander = AiZynthExpander()
    smi = "CCCCOc1ccc(CC(=O)N(C)O)cc1"
    mock_expansion_policy(TreeMolecule(parent=None, smiles=smi))
    expander.filter_policy.load("dummy.hdf5", "policy1")
    expander.filter_policy.select("policy1")

    reactions = expander.do_expansion(smi)

    assert len(reactions) == 2
    assert reactions[0][0].metadata["feasibility"] == 0.2
    assert reactions[1][0].metadata["feasibility"] == 0.2
