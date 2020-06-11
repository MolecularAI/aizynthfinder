""" Module containing helper routines for using Keras and Tensorflow
"""
import functools

from tensorflow.keras.metrics import top_k_categorical_accuracy
from tensorflow.keras.models import load_model

top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
top10_acc.__name__ = "top10_acc"

top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
top50_acc.__name__ = "top50_acc"

CUSTOM_OBJECTS = {"top10_acc": top10_acc, "top50_acc": top50_acc}


class LocalKerasModel:
    """
    A keras policy model that is executed locally.

    The size of the input vector can be determined with the len() method.

    :ivar model: the compiled model
    :vartype model: tensorflow.keras.models.Model

    :param filename: the path to a Keras checkpoint file
    :type filename: str
    """

    def __init__(self, filename):
        self.model = load_model(filename, custom_objects=CUSTOM_OBJECTS)
        try:
            self._model_dimensions = int(self.model.input.shape[1])
        except AttributeError:
            self._model_dimensions = int(self.model.input[0].shape[1])

    def __len__(self):
        return self._model_dimensions

    def predict(self, input_):
        """
        Perform a forward pass of the neural network.

        :param input_: the input vector
        :type input_: numpy.ndarray
        :return: the vector of the output layer
        :rtype: numpy.ndarray
        """
        return self.model.predict(input_)
