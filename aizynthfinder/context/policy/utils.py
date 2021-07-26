""" Module containing helper routines for policies
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np


if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Union, Any
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction


def _make_fingerprint(
    obj: Union[TreeMolecule, RetroReaction], model: Any
) -> np.ndarray:
    fingerprint = obj.fingerprint(radius=2, nbits=len(model))
    return fingerprint.reshape([1, len(model)])
