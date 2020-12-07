""" Module containing helper routines for using Keras and Tensorflow models
"""
import functools
import numpy as np
import os

import requests
import grpc
import tensorflow as tf
from google.protobuf.json_format import MessageToDict
from tensorflow_serving.apis import (
    predict_pb2,
    get_model_metadata_pb2,
    prediction_service_pb2_grpc,
)
from tensorflow.keras.metrics import top_k_categorical_accuracy
from tensorflow.keras.models import load_model as load_keras_model

from aizynthfinder.utils.logging import logger


top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
top10_acc.__name__ = "top10_acc"

top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
top50_acc.__name__ = "top50_acc"

CUSTOM_OBJECTS = {"top10_acc": top10_acc, "top50_acc": top50_acc, "tf": tf}

_logger = logger()

TF_SERVING_HOST = os.environ.get("TF_SERVING_HOST")
TF_SERVING_REST_PORT = os.environ.get("TF_SERVING_REST_PORT")
TF_SERVING_GRPC_PORT = os.environ.get("TF_SERVING_GRPC_PORT")


def load_model(source, key):
    """
    Load model from a configuration specification.

    Tries to load:
      1. A Tensorflow server through gRPC
      2. A Tensorflow server through REST API
      3. A local model

    :param source: if fallbacks to a local model, this is the filename
    :type source: str
    :param key: when connecting to Tensrflow server this is the model name
    :type key: str
    :return: a model object with a predict object
    :rtype: str
    """
    try:
        return ExternalModelViaGRPC(key)
    except ExternalModelAPIError:
        pass
    try:
        return ExternalModelViaREST(key)
    except ExternalModelAPIError:
        pass
    return LocalKerasModel(source)


class LocalKerasModel:
    """
    A keras policy model that is executed locally.

    The size of the input vector can be determined with the len() method.

    :ivar model: the compiled model
    :vartype model: tensorflow.keras.models.Model
    :ivar output_size: the length of the output vector
    :vartype output_size: int

    :param filename: the path to a Keras checkpoint file
    :type filename: str
    """

    def __init__(self, filename):
        self.model = load_keras_model(filename, custom_objects=CUSTOM_OBJECTS)
        try:
            self._model_dimensions = int(self.model.input.shape[1])
        except AttributeError:
            self._model_dimensions = int(self.model.input[0].shape[1])
        self.output_size = int(self.model.output.shape[1])

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


class ExternalModelAPIError(Exception):
    """
    Custom error type to signal failure in External model.
    """


def log_and_reraise_exceptions(method):
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        try:
            return method(*args, **kwargs)
        except Exception as e:
            msg = "Error when requesting aizynthfinder model API"
            _logger.error("%s: %s", msg, e)
            raise ExternalModelAPIError(msg)

    return wrapper


class ExternalModelViaREST:
    """
    A policy implementation using TF Serving via REST API.

    :ivar model_url: model's url
    :vartype model_url: str
    :ivar sig_def: model's signature definition (input/output)
    :vartype sig_def: dict

    :param name: the name of model (must corespond to config file)
    :type name: str
    """

    def __init__(self, name):
        self.model_url = self._get_model_url(name)
        self.sig_def = self._get_sig_def()

    def _get_model_url(self, name):
        warning = f"Failed to get predict url for external model {name}"
        if not TF_SERVING_HOST:
            _logger.warning(warning)
            raise ExternalModelAPIError("Host not set for model {name}")
        if not TF_SERVING_REST_PORT:
            _logger.warning(warning)
            raise ExternalModelAPIError("REST port not set for model {name}")
        return f"http://{TF_SERVING_HOST}:{TF_SERVING_REST_PORT}/v1/models/{name}"

    def _get_sig_def(self):
        res = self._handle_rest_api_request("GET", self.model_url + "/metadata")
        return res["metadata"]["signature_def"]["signature_def"]["serving_default"]

    def _handle_rest_api_request(self, method, url, *args, **kwargs):
        try:
            res = requests.request(method, url, *args, **kwargs)
        except Exception as e:
            msg = "Error when requesting aizynthfinder model API"
            _logger.error("%s: %s", msg, e)
            raise ExternalModelAPIError(msg)
        else:
            if res.status_code != 200 or (
                res.headers["Content-Type"] != "application/json"
            ):
                msg = "Unexpected response from aizynthfinder model API"
                _logger.error("%s: %s\n%s", msg, res.status_code, res.text)
                raise ExternalModelAPIError(msg)
            return res.json()

    def predict(self, inputs):
        """
        Get prediction from model.

        :param inputs: the input vector or list of vectors
        :type inputs: numpy.ndarray or list of numpy.ndarray
        :return: the vector of the output layer
        :rtype: numpy.ndarray
        """
        url = self.model_url + ":predict"
        res = self._handle_rest_api_request(
            "POST", url, json=self._make_payload(inputs)
        )
        return res["outputs"]

    def _make_payload(self, inputs):
        if isinstance(inputs, np.ndarray):
            inputs = [inputs]
        data = {}
        for name, fp in zip(self.sig_def["inputs"].keys(), inputs):
            data[name] = [fp.tolist()]
        return {"inputs": data}

    def __len__(self):
        first_input_name = list(self.sig_def["inputs"].keys())[0]
        return int(
            self.sig_def["inputs"][first_input_name]["tensor_shape"]["dim"][1]["size"]
        )


class ExternalModelViaGRPC:
    """
    A policy implementation using TF Serving via REST API.

    :ivar server: model's server host
    :vartype server: str
    :ivar model_name: model's name
    :vartype model_name: str
    :ivar sig_def: model's signature definition (input/output)
    :vartype sig_def: dict

    :param name: the name of model (must corespond to config file)
    :type name: str
    """

    def __init__(self, name):
        self.server = self._get_server(name)
        self.model_name = name
        self.sig_def = self._get_sig_def()

    def _get_server(self, name):
        warning = f"Failed to get server for external model {name}"
        if not TF_SERVING_HOST:
            _logger.warning(warning)
            raise ExternalModelAPIError(f"Host not set for model {name}")
        if not TF_SERVING_GRPC_PORT:
            _logger.warning(warning)
            raise ExternalModelAPIError(f"GRPC port not set for model {name}")
        return f"{TF_SERVING_HOST}:{TF_SERVING_GRPC_PORT}"

    @log_and_reraise_exceptions
    def predict(self, inputs):
        """
        Get prediction from model.

        :param inputs: the input vector or list of vectors
        :type inputs: numpy.ndarray or list of numpy.ndarray
        :return: the vector of the output layer
        :rtype: numpy.ndarray
        """
        input_tensors = self._prepare_input_tensors(inputs)
        channel = grpc.insecure_channel(self.server)
        service = prediction_service_pb2_grpc.PredictionServiceStub(channel)
        request = predict_pb2.PredictRequest()
        request.model_spec.name = self.model_name
        for name, tensor in input_tensors.items():
            request.inputs[name].CopyFrom(tensor)
        key = list(self.sig_def["outputs"].keys())[0]
        return tf.make_ndarray(service.Predict(request, 10.0).outputs[key])

    @log_and_reraise_exceptions
    def _get_sig_def(self):
        channel = grpc.insecure_channel(self.server)
        service = prediction_service_pb2_grpc.PredictionServiceStub(channel)
        request = get_model_metadata_pb2.GetModelMetadataRequest()
        request.model_spec.name = self.model_name
        request.metadata_field.append("signature_def")
        result = MessageToDict(service.GetModelMetadata(request, 10.0))
        # close the channel so that it won't be reused after fork and fail
        channel.close()
        return result["metadata"]["signature_def"]["signatureDef"]["serving_default"]

    def _prepare_input_tensors(self, inputs):
        if isinstance(inputs, np.ndarray):
            inputs = [inputs]
        tensors = {}
        for name, fp in zip(self.sig_def["inputs"].keys(), inputs):
            size = int(self.sig_def["inputs"][name]["tensorShape"]["dim"][1]["size"])
            tensors[name] = tf.make_tensor_proto(fp, dtype=np.float32, shape=(1, size))
        return tensors

    def __len__(self):
        first_input_name = list(self.sig_def["inputs"].keys())[0]
        return int(
            self.sig_def["inputs"][first_input_name]["tensorShape"]["dim"][1]["size"]
        )
