from aizynthfinder.chem import (
    FixedRetroReaction,
    SmilesBasedRetroReaction,
    TemplatedRetroReaction,
    TreeMolecule,
    UniqueMolecule,
    hash_reactions,
)


def test_retro_reaction(get_action):
    reaction = get_action(applicable=True)

    products1 = reaction.reactants
    assert products1[0][0].smiles == "CCCCOc1ccc(CC(=O)Cl)cc1"
    assert products1[0][1].smiles == "CNO"
    assert (
        reaction.reaction_smiles()
        == "CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"
    )
    assert (
        reaction.mapped_reaction_smiles().split(">>")[0]
        == "[CH3:1][CH2:2][CH2:3][CH2:4][O:5][c:6]1[cH:7][cH:8][c:9]([CH2:10]"
        "[C:11](=[O:12])[N:13]([CH3:14])[OH:15])[cH:16][cH:17]1"
    )
    assert (
        reaction.mapped_reaction_smiles().split(">>")[1]
        == "[CH3:1][CH2:2][CH2:3][CH2:4][O:5][c:6]1[cH:7][cH:8][c:9]([CH2:10]"
        "[C:11](=[O:12])[Cl:18])[cH:16][cH:17]1.[NH:13]([CH3:14])[OH:15]"
    )

    reaction = get_action(applicable=False)
    assert not reaction.reactants


def test_reaction_failure_rdchiral(get_action, mocker):
    patched_rchiral_run = mocker.patch("aizynthfinder.chem.reaction.rdc.rdchiralRun")
    patched_rchiral_run.side_effect = RuntimeError("Oh no!")
    reaction = get_action(applicable=True)

    products = reaction.reactants
    assert not products


def test_retro_reaction_with_rdkit(get_action):
    reaction = get_action(applicable=True, use_rdchiral=False)

    products1 = reaction.reactants
    assert products1[0][0].smiles == "CCCCOc1ccc(CC(=O)Cl)cc1"
    assert products1[0][1].smiles == "CNO"
    assert (
        reaction.reaction_smiles()
        == "CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"
    )

    reaction = get_action(applicable=False)
    assert not reaction.reactants


def test_retro_reaction_fingerprint(get_action):
    reaction = get_action()

    fingerprint = reaction.fingerprint(2, 10)

    assert list(fingerprint) == [0, -1, 0, -1, -1, 0, -1, -1, 0, 0]


def test_retro_reaction_copy(get_action):
    reaction = get_action()
    _ = reaction.reactants

    copy_ = reaction.copy()

    assert isinstance(copy_, TemplatedRetroReaction)
    assert copy_.mol.smiles == reaction.mol.smiles
    assert len(copy_.reactants[0]) == len(reaction.reactants[0])
    assert copy_.index == reaction.index

    copy_ = reaction.copy(index=2)

    assert copy_.mol.smiles == reaction.mol.smiles
    assert len(copy_.reactants[0]) == len(reaction.reactants[0])
    assert copy_.index != reaction.index


def test_smiles_based_retroreaction():
    mol = TreeMolecule(smiles="CNC(C)=O", parent=None)
    reaction = SmilesBasedRetroReaction(mol, reactants_str="CC(=O)O.CN")

    assert len(reaction.reactants) == 1
    assert reaction.reactants[0][0].smiles == "CC(=O)O"
    assert reaction.reactants[0][1].smiles == "CN"
    assert reaction.reaction_smiles() == "CNC(C)=O>>CC(=O)O.CN"
    assert (
        reaction.mapped_reaction_smiles()
        == "[CH3:1][NH:2][C:3]([CH3:4])=[O:5]>>[CH3:6][C:7](=[O:8])[OH:9].[CH3:10][NH2:11]"
    )


def test_smiles_based_retroreaction_w_mapping():
    mol = TreeMolecule(smiles="CNC(C)=O", parent=None, sanitize=True)
    reaction = SmilesBasedRetroReaction(
        mol, reactants_str="[CH3:5]C(=O)O.CN", mapped_prod_smiles="[CH3:5]NC(C)=O"
    )

    assert len(reaction.reactants) == 1
    assert reaction.reactants[0][0].smiles == "CC(=O)O"
    assert reaction.reactants[0][1].smiles == "CN"
    assert reaction.reaction_smiles() == "CNC(C)=O>>CC(=O)O.CN"
    assert (
        reaction.mapped_reaction_smiles()
        == "[CH3:1][NH:2][C:3]([CH3:4])=[O:5]>>[CH3:1][C:6](=[O:7])[OH:8].[CH3:9][NH2:10]"
    )


def test_create_fixed_reaction():
    smiles = "[C:1](=[O:2])([cH3:3])[N:4][cH3:5]>>Cl[C:1](=[O:2])[cH3:3].[N:4][cH3:5]"
    mol = UniqueMolecule(smiles="N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F")

    rxn = FixedRetroReaction(mol, smiles=smiles)

    assert rxn.smiles == smiles


def test_set_reactants_list_of_list():
    mol = UniqueMolecule(smiles="N#Cc1cccc(NC(=O)c2ccc(F)cc2)c1F")
    reactant1 = UniqueMolecule(smiles="N#Cc1cccc(N)c1F")
    reactant2 = UniqueMolecule(smiles="O=C(Cl)c1ccc(F)cc1")
    rxn = FixedRetroReaction(mol)

    rxn.reactants = ((reactant1, reactant2),)

    assert rxn.reactants == ((reactant1, reactant2),)


def test_fixed_retroreaction_to_smiles_based_retroreaction():
    smiles = (
        "[CH3:4][NH:3][C:2]([CH3:5])=[O:1]>>[CH3:6][NH2:5].[CH3:4][C:2]([OH:1])=[O:3]"
    )
    mol = UniqueMolecule(smiles="CNC(C)=O")
    reactant1 = UniqueMolecule(smiles="CN")
    reactant2 = UniqueMolecule(smiles="CC(=O)O")

    rxn = FixedRetroReaction(mol, metadata={"mapped_reaction_smiles": smiles})
    rxn.reactants = ((reactant1, reactant2),)
    smiles_based_retroreaction = rxn.to_smiles_based_retroreaction()

    assert isinstance(smiles_based_retroreaction, SmilesBasedRetroReaction)
    assert smiles_based_retroreaction.mol.smiles == "[CH3:4][NH:3][C:2]([CH3:5])=[O:1]"
    assert smiles_based_retroreaction.reactants[0][0].smiles == "CN"
    assert smiles_based_retroreaction.reactants[0][1].smiles == "CC(=O)O"
    assert (
        smiles_based_retroreaction.mapped_reaction_smiles()
        == "[O:1]=[C:2]([NH:3][CH3:4])[CH3:5]>>[NH2:5][CH3:6].[OH:1][C:2](=[O:3])[CH3:4]"
    )


def test_fixed_retroreaction_to_smiles_based_retroreaction_no_metadata():
    mol = UniqueMolecule(smiles="CNC(C)=O")
    reactant1 = UniqueMolecule(smiles="CN")
    reactant2 = UniqueMolecule(smiles="CC(=O)O")

    rxn = FixedRetroReaction(mol)
    rxn.reactants = ((reactant1, reactant2),)
    smiles_based_retroreaction = rxn.to_smiles_based_retroreaction()

    assert isinstance(smiles_based_retroreaction, SmilesBasedRetroReaction)
    assert smiles_based_retroreaction.mol.smiles == "CNC(C)=O"
    assert smiles_based_retroreaction.reactants[0][0].smiles == "CN"
    assert smiles_based_retroreaction.reactants[0][1].smiles == "CC(=O)O"
    assert (
        smiles_based_retroreaction.mapped_reaction_smiles()
        == "[CH3:1][NH:2][C:3]([CH3:4])=[O:5]>>[CH3:6][NH2:7].[CH3:8][C:9](=[O:10])[OH:11]"
    )


def test_reaction_hash(setup_linear_reaction_tree):
    reaction_tree = setup_linear_reaction_tree()
    reactions = list(reaction_tree.reactions())[:4]

    hash_ = hash_reactions(reactions)

    assert hash_ == "4e4ca9d7d2fc47ed3fa43a1dfb9abcd14c58f56ff73942dd1d0b8176"

    hash_ = hash_reactions(reactions, sort=False)

    assert hash_ == "567c23da4673b8b2519aeafda9b26ae949ad3e24f570968ee5f80878"


def test_mapped_atom_bonds():
    mol = TreeMolecule(smiles="[CH3:1][NH:2][C:3](C)=[O:4]", parent=None)
    reaction = SmilesBasedRetroReaction(
        mol,
        mapped_prod_smiles="[CH3:1][NH:2][C:3](C)=[O:4]",
        reactants_str="C[C:3](=[O:4])O.[CH3:1][NH:2]",
    )

    assert reaction.reactants[0][0].mapped_atom_bonds == [(6, 3), (3, 4), (3, 7)]
    assert reaction.reactants[0][1].mapped_atom_bonds == [(1, 2)]


def test_mapped_atom_bonds_rdchiral():
    smi = "CC(C)C[C:1](=O)[N:2](C)O"
    mol = TreeMolecule(smiles=smi, parent=None)
    smarts = (
        "([#8:4]-[N;H0;D3;+0:5](-[C;D1;H3:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])"
        ">>(Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([#8:4]-[NH;D2;+0:5]-[C;D1;H3:6])"
    )
    reaction = TemplatedRetroReaction(mol, smarts=smarts)

    assert reaction.reactants[0][0].mapped_atom_bonds == [
        (1, 6),
        (6, 4),
        (4, 3),
        (4, 5),
        (1, 7),
        (1, 10),
    ]
    assert reaction.reactants[0][1].mapped_atom_bonds == [(2, 8), (2, 9)]


def test_mapped_atom_bonds_rdkit():
    smi = "CC(C)C[C:1](=O)[N:2](C)O"
    mol = TreeMolecule(smiles=smi, parent=None)
    smarts = (
        "([#8:4]-[N;H0;D3;+0:5](-[C;D1;H3:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3])"
        ">>(Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]).([#8:4]-[NH;D2;+0:5]-[C;D1;H3:6])"
    )
    reaction = TemplatedRetroReaction(mol, smarts=smarts, use_rdchiral=False)

    assert reaction.reactants[0][0].mapped_atom_bonds == [
        (10, 1),
        (1, 6),
        (1, 7),
        (4, 6),
        (3, 4),
        (4, 5),
    ]
    assert reaction.reactants[0][1].mapped_atom_bonds == [(9, 2), (2, 8)]


def test_inherit_atom_tracking_rdchiral_growing():
    mol = TreeMolecule(parent=None, smiles="CN1C(=O)C2CNCCN2c2ccccc21", sanitize=True)
    smarts = (
        "[C:2]-[NH;D2;+0:1]-[C:3]>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]"
    )
    rxn = TemplatedRetroReaction(mol, smarts=smarts)

    assert all(
        atom.GetAtomMapNum() > 0 and atom.GetAtomMapNum() < 900
        for atom in rxn.reactants[0][0].mapped_mol.GetAtoms()
    )


def test_inherit_atom_tracking_rdkit_growing():
    mol = TreeMolecule(parent=None, smiles="CN1C(=O)C2CNCCN2c2ccccc21", sanitize=True)
    smarts = (
        "[C:2]-[NH;D2;+0:1]-[C:3]>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]"
    )
    rxn = TemplatedRetroReaction(mol, smarts=smarts, use_rdchiral=False)

    assert all(
        atom.GetAtomMapNum() > 0 for atom in rxn.reactants[0][0].mapped_mol.GetAtoms()
    )


def test_inherit_atom_tracking_smiles_growing():
    mol = TreeMolecule(parent=None, smiles="c1ccc(-c2ccccc2)cc1", sanitize=True)
    rxn = SmilesBasedRetroReaction(
        mol=mol,
        reactants_str="Cl[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1.Cl[c:1]1[cH:2][cH:3][cH:4][cH:11][cH:12]1",
        mapped_prod_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][cH:8][cH:9][cH:10]2)[cH:11][cH:12]1",
    )

    assert (
        rxn.reactants[0][0].mapped_smiles
        == "[c:5]1([Cl:13])[cH:6][cH:7][cH:8][cH:9][cH:10]1"
    )
    assert (
        rxn.reactants[0][1].mapped_smiles
        == "[c:1]1([Cl:14])[cH:2][cH:3][cH:4][cH:11][cH:12]1"
    )
