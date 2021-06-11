""" Module containing all types and type imports
"""
# pylint: disable=unused-import
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

from PIL.Image import Image
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

StrDict = Dict[str, Any]
RdMol = Chem.rdchem.Mol
RdReaction = Chem.rdChemReactions.ChemicalReaction
BitVector = ExplicitBitVect
PilImage = Image
PilColor = Union[str, Tuple[int, int, int]]
FrameColors = Optional[Dict[bool, PilColor]]
