import pytest
import numpy as np
import tensorflow as tf

import aizynthfinder.utils.models as models
from aizynthfinder.utils.models import ExternalModelViaREST, ExternalModelViaGRPC


@pytest.fixture()
def setup_rest_mock(mocker):
    models.TF_SERVING_HOST = "localhost"
    models.TF_SERVING_REST_PORT = "255"
    mocked_request = mocker.patch("aizynthfinder.utils.models.requests.request")
    mocked_request.return_value.status_code = 200
    mocked_request.return_value.headers = {"Content-Type": "application/json"}

    def wrapper(response):
        if isinstance(response, list):
            mocked_request.return_value.json.side_effect = response
        else:
            mocked_request.return_value.json.return_value = response
        return mocked_request

    yield wrapper

    models.TF_SERVING_HOST = None
    models.TF_SERVING_REST_PORT = None


@pytest.fixture()
def setup_grpc_mock(mocker, signature_grpc):
    models.TF_SERVING_HOST = "localhost"
    models.TF_SERVING_GRPC_PORT = "255"
    mocker.patch("aizynthfinder.utils.models.grpc.insecure_channel")
    mocked_pred_service = mocker.patch(
        "aizynthfinder.utils.models.prediction_service_pb2_grpc.PredictionServiceStub"
    )
    mocker.patch(
        "aizynthfinder.utils.models.get_model_metadata_pb2.GetModelMetadataRequest"
    )
    mocker.patch("aizynthfinder.utils.models.predict_pb2.PredictRequest")
    mocked_message = mocker.patch("aizynthfinder.utils.models.MessageToDict")
    mocked_message.return_value = signature_grpc

    def wrapper(response=None):
        if not response:
            return
        mocked_pred_service.return_value.Predict.return_value.outputs = response

    yield wrapper

    models.TF_SERVING_HOST = None
    models.TF_SERVING_GRPC_PORT = None


@pytest.fixture()
def signature_rest():
    return {
        "metadata": {
            "signature_def": {
                "signature_def": {
                    "serving_default": {
                        "inputs": {
                            "first_layer": {
                                "tensor_shape": {"dim": [{"size": 1}, {"size": 2048}]}
                            }
                        }
                    }
                }
            }
        }
    }


@pytest.fixture()
def signature_grpc():
    return {
        "metadata": {
            "signature_def": {
                "signatureDef": {
                    "serving_default": {
                        "inputs": {
                            "first_layer": {
                                "tensorShape": {"dim": [{"size": 1}, {"size": 2048}]}
                            }
                        },
                        "outputs": {"output": None},
                    }
                }
            }
        }
    }


def test_setup_tf_rest_model(signature_rest, setup_rest_mock):
    setup_rest_mock(signature_rest)

    model = ExternalModelViaREST("dummy")

    assert len(model) == 2048


def test_predict_tf_rest_model(signature_rest, setup_rest_mock):
    responses = [signature_rest, {"outputs": [0.0, 1.0]}]
    setup_rest_mock(responses)
    model = ExternalModelViaREST("dummy")

    out = model.predict(np.zeros([1, len(model)]))

    assert list(out) == [0.0, 1.0]


def test_setup_tf_grpc_model(setup_grpc_mock):
    setup_grpc_mock()

    model = ExternalModelViaGRPC("dummy")

    assert len(model) == 2048


def test_predict_tf_grpc_model(setup_grpc_mock):
    setup_grpc_mock({"output": tf.make_tensor_proto(tf.constant([0.0, 1.0]))})
    model = ExternalModelViaGRPC("dummy")

    out = model.predict(np.zeros([1, len(model)]))

    assert list(out) == [0.0, 1.0]
