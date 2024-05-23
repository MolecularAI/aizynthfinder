""" Module containing diverse math functions, including neural network-related functions. """

import numpy as np
from aizynthfinder.utils.type_utils import Callable


# pylint: disable=invalid-name
def dense_layer_forward_pass(
    x: np.ndarray, weights: np.ndarray, bias: np.ndarray, activation: Callable
) -> np.ndarray:
    """
    Forward pass through a dense neural network layer.
    :param x: layer input
    :param weights: layer weights
    :param bias: layer bias
    :param activation: layer activation function
    :return: the layer output
    """
    x = np.matmul(x, weights) + bias
    return activation(x)


# pylint: disable=invalid-name
def rectified_linear_unit(x: np.ndarray) -> np.ndarray:
    """ReLU activation function"""
    return x * (x > 0)


# pylint: disable=invalid-name
def sigmoid(x: np.ndarray) -> np.ndarray:
    """Sigmoid activation function"""
    return 1 / (1 + np.exp(-x))


# pylint: disable=invalid-name
def softmax(x: np.ndarray) -> np.ndarray:
    """Compute softmax values for each sets of scores in x."""
    return np.exp(x) / np.sum(np.exp(x), axis=0)
