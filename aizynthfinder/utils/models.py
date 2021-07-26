""" Module containing helper routines for using Keras and Tensorflow models
"""
from __future__ import annotations
import functools
import os
from typing import TYPE_CHECKING

import numpy as np
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
from aizynthfinder.utils.exceptions import ExternalModelAPIError

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import Any, Union, Callable, List

    _ModelInput = Union[np.ndarray, List[np.ndarray]]


top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
top10_acc.__name__ = "top10_acc"  # type: ignore

top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
top50_acc.__name__ = "top50_acc"  # type: ignore

CUSTOM_OBJECTS = {"top10_acc": top10_acc, "top50_acc": top50_acc, "tf": tf}

_logger = logger()

TF_SERVING_HOST = os.environ.get("TF_SERVING_HOST")
TF_SERVING_REST_PORT = os.environ.get("TF_SERVING_REST_PORT")
TF_SERVING_GRPC_PORT = os.environ.get("TF_SERVING_GRPC_PORT")


def load_model(
    source: str, key: str, use_remote_models: bool
) -> Union["LocalKerasModel", "ExternalModelViaGRPC", "ExternalModelViaREST"]:
    """
    Load model from a configuration specification.

    If `use_remote_models` is True, tries to load:
      1. A Tensorflow server through gRPC
      2. A Tensorflow server through REST API
      3. A local model
    otherwise it just loads the local model

    :param source: if fallbacks to a local model, this is the filename
    :param key: when connecting to Tensrflow server this is the model name
    :param use_remote_models: if True will try to connect to remote model server
    :return: a model object with a predict object
    """
    if not use_remote_models:
        return LocalKerasModel(source)

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
    :ivar output_size: the length of the output vector

    :param filename: the path to a Keras checkpoint file
    """

    def __init__(self, filename: str) -> None:
        self.model = load_keras_model(filename, custom_objects=CUSTOM_OBJECTS)
        try:
            self._model_dimensions = int(self.model.input.shape[1])
        except AttributeError:
            self._model_dimensions = int(self.model.input[0].shape[1])
        self.output_size = int(self.model.output.shape[1])

    def __len__(self) -> int:
        return self._model_dimensions

    def predict(self, input_: _ModelInput) -> np.ndarray:
        """
        Perform a forward pass of the neural network.

        :param input_: the input vector
        :return: the vector of the output layer
        """
        return self.model.predict(input_)


def _log_and_reraise_exceptions(method: Callable) -> Callable:
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        try:
            return method(*args, **kwargs)
        except Exception as err:
            msg = "Error when requesting from tensorflow model API"
            _logger.error("%s: %s", msg, err)
            raise ExternalModelAPIError(msg)

    return wrapper


class ExternalModelViaREST:
    """
    A neural network model implementation using TF Serving via REST API.

    :param name: the name of model
    """

    def __init__(self, name: str) -> None:
        self._model_url = self._get_model_url(name)
        self._sig_def = self._get_sig_def()

    def __len__(self) -> int:
        first_input_name = list(self._sig_def["inputs"].keys())[0]
        return int(
            self._sig_def["inputs"][first_input_name]["tensor_shape"]["dim"][1]["size"]
        )

    def predict(self, inputs: _ModelInput) -> np.ndarray:
        """
        Get prediction from model.

        :param inputs: the input vector or list of vectors
        :return: the vector of the output layer
        """
        url = self._model_url + ":predict"
        res = self._handle_rest_api_request(
            "POST", url, json=self._make_payload(inputs)
        )
        return np.asarray(res["outputs"])

    def _get_sig_def(self) -> dict:
        res = self._handle_rest_api_request("GET", self._model_url + "/metadata")
        return res["metadata"]["signature_def"]["signature_def"]["serving_default"]

    # pylint: disable=no-self-use
    @_log_and_reraise_exceptions
    def _handle_rest_api_request(
        self, method: str, url: str, *args: Any, **kwargs: Any
    ) -> dict:
        res = requests.request(method, url, *args, **kwargs)
        if res.status_code != 200 or (
            res.headers["Content-Type"] != "application/json"
        ):
            raise ExternalModelAPIError(
                f"Unexpected response from REST API: {res.status_code}\n{res.text}"
            )

        return res.json()

    def _make_payload(self, inputs: _ModelInput) -> dict:
        if isinstance(inputs, np.ndarray):
            inputs = [inputs]
        data = {
            name: fp.tolist()
            for name, fp in zip(self._sig_def["inputs"].keys(), inputs)
        }
        return {"inputs": data}

    @staticmethod
    def _get_model_url(name: str) -> str:
        warning = f"Failed to get url of REST service for external model {name}"
        if not TF_SERVING_HOST:
            _logger.warning(warning)
            raise ExternalModelAPIError("Host not set for model {name}")
        if not TF_SERVING_REST_PORT:
            _logger.warning(warning)
            raise ExternalModelAPIError("REST port not set for model {name}")
        return f"http://{TF_SERVING_HOST}:{TF_SERVING_REST_PORT}/v1/models/{name}"


class ExternalModelViaGRPC:
    """
    A neural network model implementation using TF Serving via gRPC.

    :param name: the name of model
    """

    def __init__(self, name: str) -> None:
        self._server = self._get_server(name)
        self._model_name = name
        self._sig_def = self._get_sig_def()

    def __len__(self) -> int:
        first_input_name = list(self._sig_def["inputs"].keys())[0]
        return int(
            self._sig_def["inputs"][first_input_name]["tensorShape"]["dim"][1]["size"]
        )

    @_log_and_reraise_exceptions
    def predict(self, inputs: _ModelInput) -> np.ndarray:
        """
        Get prediction from model.

        :param inputs: the input vector or list of vectors
        :return: the vector of the output layer
        """
        input_tensors = self._make_payload(inputs)
        channel = grpc.insecure_channel(self._server)
        service = prediction_service_pb2_grpc.PredictionServiceStub(channel)
        request = predict_pb2.PredictRequest()
        request.model_spec.name = self._model_name
        for name, tensor in input_tensors.items():
            request.inputs[name].CopyFrom(tensor)
        key = list(self._sig_def["outputs"].keys())[0]
        return tf.make_ndarray(service.Predict(request, 10.0).outputs[key])

    @_log_and_reraise_exceptions
    def _get_sig_def(self) -> dict:
        channel = grpc.insecure_channel(self._server)
        service = prediction_service_pb2_grpc.PredictionServiceStub(channel)
        request = get_model_metadata_pb2.GetModelMetadataRequest()
        request.model_spec.name = self._model_name
        request.metadata_field.append("signature_def")
        result = MessageToDict(service.GetModelMetadata(request, 10.0))
        # close the channel so that it won't be reused after fork and fail
        channel.close()
        return result["metadata"]["signature_def"]["signatureDef"]["serving_default"]

    def _make_payload(self, inputs: _ModelInput) -> dict:
        if isinstance(inputs, np.ndarray):
            inputs = [inputs]
        tensors = {}
        for name, fp_ in zip(self._sig_def["inputs"].keys(), inputs):
            size = int(self._sig_def["inputs"][name]["tensorShape"]["dim"][1]["size"])
            tensors[name] = tf.make_tensor_proto(fp_, dtype=np.float32, shape=(1, size))
        return tensors

    @staticmethod
    def _get_server(name: str) -> str:
        warning = f"Failed to get gRPC server for external model {name}"
        if not TF_SERVING_HOST:
            _logger.warning(warning)
            raise ExternalModelAPIError(f"Host not set for model {name}")
        if not TF_SERVING_GRPC_PORT:
            _logger.warning(warning)
            raise ExternalModelAPIError(f"GRPC port not set for model {name}")
        return f"{TF_SERVING_HOST}:{TF_SERVING_GRPC_PORT}"
