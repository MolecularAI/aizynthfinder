""" Module containing classes to deal with Molecules and Reactions - mostly wrappers around rdkit routines.
"""
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdchiral import main as rdc


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
    def inchi_key(self):
        """
        The inchi key representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi key
        :rtype: str
        """
        if not self._inchi_key:
            self.sanitize()
            self._inchi_key = Chem.MolToInchiKey(self.rd_mol)
        return self._inchi_key

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

    def make_unique(self):
        """
        Returns an instance of the UniqueMolecule class that
        is representing the same molecule but is not hashable or comparable.

        :return: the unique molecule
        :rtype: UniqueMolecule
        """
        return UniqueMolecule(rd_mol=self.rd_mol)

    def sanitize(self):
        """
        Sanitizes the molecule if it has not been done before.

        :raises MoleculeException: if the molecule could not be sanitized
        """
        if self._is_sanitized:
            return

        try:
            AllChem.SanitizeMol(self.rd_mol)
        except:  # noqa, there could be many reasons why the molecule cannot be sanitized
            raise MoleculeException("Unable to sanitize molecule")

        self.smiles = Chem.MolToSmiles(self.rd_mol)
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

    :ivar mol: the TreeMolecule object that this reaction is applied to
    :vartype mol: TreeMolecule
    :ivar smarts: the SMARTS representation of the reaction
    :vartype smarts: str
    :ivar index: a unique index of this reaction,
                 to count for the fact that a reaction can produce more than one outcome
    :vartype index: int
    :ivar metadata: meta data associated with the reaction
    :vartype metadata: dict

    :param mol: the molecule
    :type mol: TreeMolecule
    :param smarts: the SMARTS fragment
    :type smarts: str
    :param index: the index, defaults to 0
    :type index: int, optional
    :ivar metadata: some meta data
    :vartype metadata: dict, optional
    """

    def __init__(self, mol, smarts, index=0, metadata=None):
        self.mol = mol
        self.smarts = smarts
        self.index = index
        self.metadata = metadata or dict()
        self._reactants = None
        self._rd_reaction = None
        self._smiles = None

    def __str__(self):
        return f"template {self.smarts} on molecule {self.mol.smiles}"

    @property
    def reactants(self):
        """
        Returns the reactant molcules.
        Apply the reaction if necessary.

        :return: the products of the reaction
        :rtype: tuple of tuple of TreeMolecule
        """
        if not self._reactants:
            self.apply()
        return self._reactants

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
        """
        Apply a reactions smarts to a molecule and return the products (reactants for retro templates)

        Will try to sanitize the reactants, and if that fails it will not return that molecule

        :return: the products of the reaction
        :rtype: tuple of tuple of TreeMolecule
        """
        reaction = rdc.rdchiralReaction(self.smarts)
        rct = rdc.rdchiralReactants(self.mol.smiles)
        reactants = rdc.rdchiralRun(reaction, rct)

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
        self._reactants = tuple(outcomes)

        return self._reactants

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
        reactants = ".".join(reactant.smiles for reactant in self.reactants[self.index])
        return "%s>>%s" % (reactants, self.mol.smiles)
