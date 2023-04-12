from typing import Dict, List

import numpy as np
import pytest
import pytest_mock
from aizynthfinder.utils import models


def test_local_onnx_model_predict(mock_onnx_model: pytest_mock.MockerFixture) -> None:
    onnx_model = models.LocalOnnxModel("test_model.onnx")
    output = onnx_model.predict(np.array([1]))
    expected_output = np.array([[0.2, 0.7, 0.1]])

    assert np.array_equal(output, expected_output)


def test_local_onnx_model_length(mock_onnx_model: pytest_mock.MockerFixture) -> None:
    onnx_model = models.LocalOnnxModel("test_model.onnx")
    output = len(onnx_model)
    expected_output = 3

    assert output == expected_output


def test_local_onnx_model_output_size(
    mock_onnx_model: pytest_mock.MockerFixture,
) -> None:
    onnx_model = models.LocalOnnxModel("test_model.onnx")
    output = onnx_model.output_size
    expected_output = 3

    assert output == expected_output
