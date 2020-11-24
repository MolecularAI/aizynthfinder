""" Module containing classes to deal with Molecules and Reactions - mostly wrappers around rdkit routines.
"""
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdchiral import main as rdc

from aizynthfinder.utils.logging import logger


class MoleculeException(Exception):
    """ An exception that is raised by the Molecule class
    """


class Molecule:
    """
    A base class for molecules. Encapsulate an RDKit mol object and
    functions that can be applied to such a molecule.

    The objects of this class is hashable by the inchi key and hence
    comparable with the equality operator.

    :ivar rd_mol: the RDkit mol object that is encapsulated
    :vartype rd_mol: rdkit.Chem.rdchem.Mol
    :ivar smiles: the SMILES representation of the molecule
    :vartype smiles: str

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :type rd_mol: rdkit.Chem.rdchem.Mol or None, optional
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :type smiles: str or None, optional
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :type sanitize: bool, optional
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(self, rd_mol=None, smiles=None, sanitize=False):
        if not rd_mol and not smiles:
            raise MoleculeException(
                "Need to provide either a rdkit Mol object or smiles to create a molecule"
            )

        if rd_mol:
            self.rd_mol = rd_mol
            self.smiles = Chem.MolToSmiles(rd_mol)
        else:
            self.smiles = smiles
            self.rd_mol = Chem.MolFromSmiles(smiles, sanitize=False)

        self._inchi_key = None
        self._inchi = None
        self._fingerprints = {}
        self._is_sanitized = False

        if sanitize:
            self.sanitize()

    def __hash__(self):
        return hash(self.inchi_key)

    def __eq__(self, other):
        return self.inchi_key == other.inchi_key

    def __str__(self):
        return self.smiles

    @property
    def inchi(self):
        """
        The inchi representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi
        :rtype: str
        """
        if not self._inchi:
            self.sanitize(raise_exception=False)
            self._inchi = Chem.MolToInchi(self.rd_mol)
        return self._inchi

    @property
    def inchi_key(self):
        """
        The inchi key representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi key
        :rtype: str
        """
        if not self._inchi_key:
            self.sanitize(raise_exception=False)
            self._inchi_key = Chem.MolToInchiKey(self.rd_mol)
        return self._inchi_key

    def basic_compare(self, other):
        """
        Compare this molecule to another but only to
        the basic part of the inchi key, thereby ignoring stereochemistry etc

        :param other: the molecule to compare to
        :type other: Molecule
        :return: True if chemical formula and connectivity is the same
        :rtype: bool
        """
        return self.inchi_key[:14] == other.inchi_key[:14]

    def fingerprint(self, radius, nbits=None):
        """
        Returns the Morgan fingerprint of the molecule

        :param radius: the radius of the fingerprint
        :type radius: int
        :param nbits: the length of the fingerprint. If not given it will use RDKit default, defaults to None
        :type nbits: int, optional
        :return: the fingerprint
        :rtype: numpy.ndarray
        """
        if nbits:
            key = (radius, nbits)
        else:
            key = (radius,)

        if key not in self._fingerprints:
            self.sanitize()
            bitvect = AllChem.GetMorganFingerprintAsBitVect(self.rd_mol, *key)
            array = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(bitvect, array)
            self._fingerprints[key] = array

        return self._fingerprints[key]

    def has_atom_mapping(self):
        """
        Determines if a the molecule has atom mappings

        :return: True if at least one atom has a mapping
        :rtype: bool
        """
        for atom in self.rd_mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                return True
        return False

    def make_unique(self):
        """
        Returns an instance of the UniqueMolecule class that
        is representing the same molecule but is not hashable or comparable.

        :return: the unique molecule
        :rtype: UniqueMolecule
        """
        return UniqueMolecule(rd_mol=self.rd_mol)

    def remove_atom_mapping(self):
        """
        Remove all mappings of the atoms and update the smiles
        """
        for atom in self.rd_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        self.smiles = Chem.MolToSmiles(self.rd_mol)

    def sanitize(self, raise_exception=True):
        """
        Sanitizes the molecule if it has not been done before.

        :param raise_exception: if True will raise exception on failed sanitaion
        :type raise_exception: bool
        :raises MoleculeException: if the molecule could not be sanitized
        """
        if self._is_sanitized:
            return

        try:
            AllChem.SanitizeMol(self.rd_mol)
        except:  # noqa, there could be many reasons why the molecule cannot be sanitized
            if raise_exception:
                raise MoleculeException(f"Unable to sanitize molecule ({self.smiles})")
            self.rd_mol = Chem.MolFromSmiles(self.smiles, sanitize=False)
            return

        self.smiles = Chem.MolToSmiles(self.rd_mol)
        self._inchi = None
        self._inchi_key = None
        self._fingerprints = {}
        self._is_sanitized = True


class TreeMolecule(Molecule):
    """
    A special molecule that keeps a reference to a parent molecule.

    If the class is intantiated without specifying the `transform` argument,
    it is computed by increasing the value of the `parent.transform` variable.

    :ivar parent: parent molecule
    :vartype parent: Molecule
    :ivar transform: a numerical number corresponding to the depth in the tree
    :vartype transform: int

    :param parent: a TreeMolecule object that is the parent
    :type parent: TreeMolecule
    :param transform: the transform value, defaults to None
    :type transform: int, optional
    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :type rd_mol: mol or None, optional
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :type smiles: str or None, optional
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :type sanitize: bool, optional
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(
        self, parent, transform=None, rd_mol=None, smiles=None, sanitize=False
    ):
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)
        self.parent = parent
        if transform is None and parent and parent.transform is not None:
            self.transform = parent.transform + 1
        else:
            self.transform = transform


class UniqueMolecule(Molecule):
    """
    A special molecule with the hash set to the id of the object.
    Therefore no two instances of this class will be comparable.

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :type rd_mol: mol or None, optional
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :type smiles: str or None, optional
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :type sanitize: bool, optional
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(self, rd_mol=None, smiles=None, sanitize=False):
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)

    def __hash__(self):
        return id(self)

    def __eq__(self, _):
        return False


class Reaction:
    """
    An abstraction of a reaction. Encapsulate an RDKit reaction object and
    functions that can be applied to such a reaction.

    :ivar mols: the Molecule objects that this reaction are applied to
    :vartype mols: list of Molecule
    :ivar smarts: the SMARTS representation of the reaction
    :vartype smarts: str
    :ivar index: a unique index of this reaction,
                 to count for the fact that a reaction can produce more than one outcome
    :vartype index: int
    :ivar metadata: meta data associated with the reaction
    :vartype metadata: dict

    :param mols: the molecules
    :type mols: list of Molecule
    :param smarts: the SMARTS fragment
    :type smarts: str
    :param index: the index, defaults to 0
    :type index: int, optional
    :ivar metadata: some meta data
    :vartype metadata: dict, optional
    """

    def __init__(self, mols, smarts, index=0, metadata=None):
        self.mols = mols
        self.smarts = smarts
        self.index = index
        self.metadata = metadata or dict()
        self._products = None
        self._rd_reaction = None
        self._smiles = None

    def __str__(self):
        return f"template {self.smarts} on molecule {self.mols[0].smiles}"

    @property
    def products(self):
        """
        Returns the product molcules.
        Apply the reaction if necessary.

        :return: the products of the reaction
        :rtype: tuple of tuple of Molecule
        """
        if not self._products:
            self.apply()
        return self._products

    @property
    def rd_reaction(self):
        """
        The reaction as a RDkit reaction object

        :return: the reaction object
        :rtype: rdkit.Chem.rdChemReactions.ChemicalReaction
        """
        if not self._rd_reaction:
            self._rd_reaction = AllChem.ReactionFromSmarts(self.smarts)
        return self._rd_reaction

    @property
    def smiles(self):
        """
        The reaction as a SMILES

        :return: the SMILES
        :rtype: str
        """
        if self._smiles is None:
            self._smiles = AllChem.ReactionToSmiles(self.rd_reaction)
        return self._smiles

    def apply(self):
        num_rectantant_templates = self.rd_reaction.GetNumReactantTemplates()
        reactants = tuple(mol.rd_mol for mol in self.mols[:num_rectantant_templates])
        products_list = self.rd_reaction.RunReactants(reactants)

        outcomes = []
        for products in products_list:
            try:
                mols = tuple(Molecule(rd_mol=mol, sanitize=True) for mol in products)
            except MoleculeException:
                pass
            else:
                outcomes.append(mols)
        self._products = tuple(outcomes)

        return self._products

    def rd_reaction_from_smiles(self):
        """
        The reaction as a RDkit reaction object but created from the reaction smiles
        instead of the SMARTS of the template.

        :return: the reaction object
        :rtype: rdkit.Chem.rdChemReactions.ChemicalReaction
        """
        return AllChem.ReactionFromSmarts(self.reaction_smiles(), useSmiles=True)

    def reaction_smiles(self):
        """
        Get the reaction SMILES, i.e. the SMILES of the reactants and products joined together

        :return: the SMILES
        :rtype: str
        """
        reactants = ".".join(mol.smiles for mol in self.mols)
        products = ".".join(mol.smiles for mol in self.products[self.index])
        return "%s>>%s" % (reactants, products)


class RetroReaction(Reaction):
    """
    A retrosynthesis reaction. Only a single molecule is the reactant.

    :ivar mol: the TreeMolecule object that this reaction is applied to
    :vartype mol: TreeMolecule

    :param mol: the molecule
    :type mol: TreeMolecule
    :param smarts: the SMARTS fragment
    :type smarts: str
    :param index: the index, defaults to 0
    :type index: int, optional
    """

    def __init__(self, mol, smarts, index=0, metadata=None):
        super().__init__([mol], smarts, index, metadata)
        self.mol = mol

    @classmethod
    def from_reaction_smiles(cls, smiles, smarts):
        """
        Construct a retro reaction by parsing a reaction smiles.

        Note that applying reaction does not necessarily give the
        same outcome.

        :param smiles: the reaction smiles
        :type smiles: str
        :param smarts: the SMARTS of the reaction
        :type smarts: str
        :return: the constructed reaction object
        :rtype: RetroReaction
        """
        mol_smiles, reactants_smiles = smiles.split(">>")
        mol = TreeMolecule(parent=None, smiles=mol_smiles)
        reaction = RetroReaction(mol, smarts=smarts)
        reaction._products = (
            [
                TreeMolecule(parent=mol, smiles=smiles)
                for smiles in reactants_smiles.split(".")
            ],
        )
        return reaction

    @property
    def reactants(self):
        """
        Returns the reactant molcules.
        Apply the reaction if necessary.

        :return: the products of the reaction
        :rtype: tuple of tuple of TreeMolecule
        """
        return self.products

    def apply(self):
        """
        Apply a reactions smarts to a molecule and return the products (reactants for retro templates)

        Will try to sanitize the reactants, and if that fails it will not return that molecule

        :return: the products of the reaction
        :rtype: tuple of tuple of TreeMolecule
        """
        reaction = rdc.rdchiralReaction(self.smarts)
        rct = rdc.rdchiralReactants(self.mol.smiles)
        try:
            reactants = rdc.rdchiralRun(reaction, rct)
        except RuntimeError as err:
            logger().debug(
                f"Runtime error in RDChiral with template {self.smarts} on {self.mol.smiles}\n{err}"
            )
            reactants = []

        # Turning rdchiral outcome into rdkit tuple of tuples to maintain compatibility
        outcomes = []
        for reactant_str in reactants:
            smiles_list = reactant_str.split(".")
            try:
                rct = tuple(
                    TreeMolecule(parent=self.mol, smiles=smi, sanitize=True)
                    for smi in smiles_list
                )
            except MoleculeException:
                pass
            else:
                outcomes.append(rct)
        self._products = tuple(outcomes)

        return self._products

    def copy(self, index=None):
        """
        Shallow copy of this instance.

        :param index: new index, defaults to None
        :type index: int, optional
        :return: the copy
        :rtype: RetroReaction
        """
        index = index if index is not None else self.index
        new_reaction = RetroReaction(self.mol, self.smarts, index, self.metadata)
        new_reaction._products = tuple(products for products in self._products)
        new_reaction._rd_reaction = self._rd_reaction
        new_reaction._smiles = self._smiles
        return new_reaction

    def fingerprint(self, radius, nbits=None):
        """
        Returns a difference fingerprint

        :param radius: the radius of the fingerprint
        :type radius: int
        :param nbits: the length of the fingerprint. If not given it will use RDKit default, defaults to None
        :type nbits: int, optional
        :return: the fingerprint
        :rtype: numpy.ndarray
        """
        product_fp = self.mol.fingerprint(radius, nbits)
        reactants_fp = sum(
            mol.fingerprint(radius, nbits) for mol in self.reactants[self.index]
        )
        return product_fp - reactants_fp
