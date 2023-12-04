""" Module containing diverse math functions, including neural network-related functions. """
import numpy as np

# pylint: disable=invalid-name
def softmax(x: np.ndarray) -> np.ndarray:
    """Compute softmax values for each sets of scores in x."""
    return np.exp(x) / np.sum(np.exp(x), axis=0)
