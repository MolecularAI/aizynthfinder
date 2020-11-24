import os

import pytest
import pandas as pd
import numpy as np
from scipy import sparse

from aizynthfinder.training.utils import (
    Config,
    split_and_save_data,
    smiles_to_fingerprint,
)
from aizynthfinder.training.keras_models import RolloutModelSequence


@pytest.fixture
def default_config():
    return Config()


@pytest.fixture
def rollout_model_sequence(mocker, default_config):
    mocked_load_npz = mocker.patch(
        "aizynthfinder.training.keras_models.sparse.load_npz"
    )
    input_array = sparse.csr_matrix(np.zeros([9000, 10]))
    label_array = sparse.csr_matrix(np.zeros([9000, 100]))
    mocked_load_npz.side_effect = [input_array, label_array]

    return RolloutModelSequence(default_config, "training")


def test_empty_config(default_config, write_yaml):
    filename = write_yaml({})

    config = Config(filename)

    assert config._config == default_config._config


def test_update_single_setting(default_config, write_yaml):
    filename = write_yaml({"fingerprint_len": 10})

    config = Config(filename)

    assert config["fingerprint_len"] == 10
    assert config["template_occurance"] == default_config["template_occurance"]


def test_update_nested_setting(default_config, write_yaml):
    filename = write_yaml(
        {"split_size": {"training": 0.8, "testing": 0.1, "validation": 0.1}}
    )

    config = Config(filename)

    assert config["template_occurance"] == default_config["template_occurance"]
    assert config["split_size"]["training"] == 0.8
    assert config["split_size"]["testing"] == 0.1
    assert config["split_size"]["validation"] == 0.1


def test_update_invalid_setting(default_config, write_yaml):
    filename = write_yaml(
        {"fingerprint_len": {"training": 0.8, "testing": 0.1, "validation": 0.1}}
    )

    config = Config(filename)

    assert config["fingerprint_len"] == default_config["fingerprint_len"]


def test_config_filename(default_config):

    filename = default_config.filename("raw_library")
    assert filename.startswith(default_config["output_path"])
    assert filename.endswith(default_config["file_postfix"]["raw_library"])

    filename = default_config.filename("something")
    assert filename.startswith(default_config["output_path"])
    assert filename.endswith("something")


def test_split_and_save_data_frame(mocker, tmpdir, default_config):
    default_config["output_path"] = str(tmpdir)
    default_config["file_prefix"] = "dummy"
    filename_train = str(
        tmpdir / "dummy" + default_config["file_postfix"]["training_library"]
    )
    filename_valid = str(
        tmpdir / "dummy" + default_config["file_postfix"]["validation_library"]
    )
    filename_test = str(
        tmpdir / "dummy" + default_config["file_postfix"]["testing_library"]
    )
    data = pd.DataFrame.from_dict({"one": np.zeros(100), "two": np.ones(100)})

    split_and_save_data(data, "library", default_config)

    assert os.path.exists(filename_train)
    assert os.path.exists(filename_valid)
    assert os.path.exists(filename_test)

    data_read = pd.read_csv(filename_train, header=None, names=["one", "two"])
    assert len(data_read) == 90

    data_read = pd.read_csv(filename_valid, header=None, names=["one", "two"])
    assert len(data_read) == 5

    data_read = pd.read_csv(filename_test, header=None, names=["one", "two"])
    assert len(data_read) == 5


def test_split_and_save_data_ndarray(mocker, tmpdir, default_config):
    default_config["output_path"] = str(tmpdir)
    default_config["file_prefix"] = "dummy"
    filename_train = str(
        tmpdir / "dummy" + default_config["file_postfix"]["training_inputs"]
    )
    filename_valid = str(
        tmpdir / "dummy" + default_config["file_postfix"]["validation_inputs"]
    )
    filename_test = str(
        tmpdir / "dummy" + default_config["file_postfix"]["testing_inputs"]
    )
    data = np.ones([100, 2])

    split_and_save_data(data, "inputs", default_config)

    assert os.path.exists(filename_train)
    assert os.path.exists(filename_valid)
    assert os.path.exists(filename_test)

    data_read = np.load(filename_train)["arr_0"]
    assert len(data_read) == 90

    data_read = np.load(filename_valid)["arr_0"]
    assert len(data_read) == 5

    data_read = np.load(filename_test)["arr_0"]
    assert len(data_read) == 5


def test_split_and_save_data_sparse(default_config, mocker, tmpdir):
    default_config["output_path"] = str(tmpdir)
    default_config["file_prefix"] = "dummy"
    filename_train = str(
        tmpdir / "dummy" + default_config["file_postfix"]["training_inputs"]
    )
    filename_valid = str(
        tmpdir / "dummy" + default_config["file_postfix"]["validation_inputs"]
    )
    filename_test = str(
        tmpdir / "dummy" + default_config["file_postfix"]["testing_inputs"]
    )
    data = sparse.csr_matrix(np.ones([100, 2]))

    split_and_save_data(data, "inputs", default_config)

    assert os.path.exists(filename_train)
    assert os.path.exists(filename_valid)
    assert os.path.exists(filename_test)

    data_read = sparse.load_npz(str(filename_train))
    assert data_read.shape[0] == 90

    data_read = sparse.load_npz(str(filename_valid))
    assert data_read.shape[0] == 5

    data_read = sparse.load_npz(str(filename_test))
    assert data_read.shape[0] == 5


def test_smiles_to_fingerprint(default_config):
    default_config["fingerprint_len"] = 10

    fingerprint = smiles_to_fingerprint(["O"], default_config)

    assert sum(fingerprint) == 1


def test_rollout_model_sequence_loading(rollout_model_sequence):

    assert rollout_model_sequence.input_dim == 10
    assert rollout_model_sequence.output_dim == 100


def test_rollout_model_sequence_slicing(rollout_model_sequence, default_config):
    seq = rollout_model_sequence

    xbatch, ybatch = seq[1]

    assert xbatch.shape[0] == default_config["batch_size"]
    assert ybatch.shape[0] == default_config["batch_size"]
    assert xbatch.shape[1] == rollout_model_sequence.input_dim
    assert ybatch.shape[1] == rollout_model_sequence.output_dim
