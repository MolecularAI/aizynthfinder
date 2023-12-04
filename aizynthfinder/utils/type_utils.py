""" Module containing all types and type imports
"""
# pylint: disable=unused-import
from typing import Callable  # noqa
from typing import Iterable  # noqa
from typing import List  # noqa
from typing import Sequence  # noqa
from typing import Set  # noqa
from typing import TypeVar  # noqa
from typing import Any, Dict, Optional, Tuple, Union

from PIL.Image import Image
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

StrDict = Dict[str, Any]
RdMol = Chem.rdchem.Mol
RdReaction = Chem.rdChemReactions.ChemicalReaction
BitVector = ExplicitBitVect
PilImage = Image
PilColor = Union[str, Tuple[int, int, int]]
FrameColors = Optional[Dict[bool, PilColor]]
