from aizynthfinder.chem import (
    Molecule,
    UniqueMolecule,
    TreeMolecule,
    Reaction,
    TemplatedRetroReaction,
    FixedRetroReaction,
    SmilesBasedRetroReaction,
    hash_reactions,
)
from aizynthfinder.reactiontree import ReactionTree


def test_fwd_reaction():
    mol1 = Molecule(smiles="CC(=O)O", sanitize=True)
    mol2 = Molecule(smiles="NC", sanitize=True)
    reaction = Reaction(
        mols=[mol1, mol2], smarts="[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
    )

    products = reaction.apply()

    assert len(products) == 1
    assert len(products[0]) == 1
    assert products[0][0].smiles == "CNC(C)=O"
    assert reaction.reaction_smiles() == "CC(=O)O.CN>>CNC(C)=O"


def test_retro_reaction(get_action):
    reaction = get_action(applicable=True)

    products1 = reaction.reactants
    assert products1[0][0].smiles == "CCCCOc1ccc(CC(=O)Cl)cc1"
    assert products1[0][1].smiles == "CNO"
    assert (
        reaction.reaction_smiles()
        == "CCCCOc1ccc(CC(=O)N(C)O)cc1>>CCCCOc1ccc(CC(=O)Cl)cc1.CNO"
    )

    reaction = get_action(applicable=False)
    assert not reaction.reactants


def test_reaction_failure_rdchiral(get_action, mocker):
    patched_rchiral_run = mocker.patch("aizynthfinder.chem.reaction.rdc.rdchiralRun")
    patched_rchiral_run.side_effect = RuntimeError("Oh no!")
    reaction = get_action(applicable=True)

    products = reaction.reactants
    assert not products


def test_retro_reaction_fingerprint(get_action):
    reaction = get_action()

    fp = reaction.fingerprint(2, 10)

    assert list(fp) == [0, -1, 0, -1, -1, 0, -1, -1, 0, 0]


def test_retro_reaction_copy(get_action):
    mol = TreeMolecule(parent=None, smiles="CCCCOc1ccc(CC(=O)N(C)O)cc1")
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
    assert reaction.smiles == "CNC(C)=O>>CC(=O)O.CN"


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


def test_reaction_hash(setup_linear_reaction_tree):
    rt = setup_linear_reaction_tree()
    reactions = list(rt.reactions())[:4]

    hash_ = hash_reactions(reactions)

    assert hash_ == "4e4ca9d7d2fc47ed3fa43a1dfb9abcd14c58f56ff73942dd1d0b8176"

    hash_ = hash_reactions(reactions, sort=False)

    assert hash_ == "567c23da4673b8b2519aeafda9b26ae949ad3e24f570968ee5f80878"
