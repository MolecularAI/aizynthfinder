from typing import (
    Any,
    Callable,  # noqa
    Dict,
    List,  # noqa
    Iterable,  # noqa
    Optional,
    Sequence,  # noqa
    Set,  # noqa
    Tuple,
    TypeVar,  # noqa
    Union,
)

from rdkit import Chem
from PIL.Image import Image
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

StrDict = Dict[str, Any]
RdMol = Chem.rdchem.Mol
RdReaction = Chem.rdChemReactions.ChemicalReaction
BitVector = ExplicitBitVect
PilImage = Image
PilColor = Union[str, Tuple[int, int, int]]
FrameColors = Optional[Dict[bool, PilColor]]
