import os

import pytest
import pandas as pd
import numpy as np
from scipy import sparse

from aizynthfinder.training.utils import (
    Config,
    create_reactants_molecules,
    reverse_template,
    reaction_hash,
    split_and_save_data,
    smiles_to_fingerprint,
    is_sanitizable,
    reaction_to_fingerprints,
)
from aizynthfinder.training.keras_models import (
    ExpansionModelSequence,
    FilterModelSequence,
)
from aizynthfinder.training.make_false_products import (
    strict_application,
    random_application,
    recommender_application,
)
from aizynthfinder.chem import Molecule


@pytest.fixture
def default_config():
    return Config()


@pytest.fixture
def filter_model_sequence(mocker, default_config):
    mocked_load_npz = mocker.patch(
        "aizynthfinder.training.keras_models.sparse.load_npz"
    )
    input_array = sparse.csr_matrix(np.zeros([9000, 10]))
    input_array2 = sparse.csr_matrix(np.zeros([9000, 10]))
    mocked_load_npz.side_effect = [input_array, ValueError(), input_array2]

    label_array = np.zeros([9000, 100])
    mocked_np_load = mocker.patch("aizynthfinder.training.keras_models.np.load")
    mocked_np_load.return_value = {"arr_0": label_array}

    return FilterModelSequence(default_config, "training")


@pytest.fixture
def expansion_model_sequence(mocker, default_config):
    mocked_load_npz = mocker.patch(
        "aizynthfinder.training.keras_models.sparse.load_npz"
    )
    input_array = sparse.csr_matrix(np.zeros([9000, 10]))
    label_array = sparse.csr_matrix(np.zeros([9000, 100]))
    mocked_load_npz.side_effect = [input_array, label_array]

    return ExpansionModelSequence(default_config, "training")


def test_empty_config(default_config, write_yaml):
    filename = write_yaml({})

    config = Config(filename)

    assert config._config == default_config._config


def test_update_single_setting(default_config, write_yaml):
    filename = write_yaml({"fingerprint_len": 10})

    config = Config(filename)

    assert config["fingerprint_len"] == 10
    assert config["template_occurrence"] == default_config["template_occurrence"]


def test_update_nested_setting(default_config, write_yaml):
    filename = write_yaml(
        {"split_size": {"training": 0.8, "testing": 0.1, "validation": 0.1}}
    )

    config = Config(filename)

    assert config["template_occurrence"] == default_config["template_occurrence"]
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


def test_split_and_save_data_frame(tmpdir, default_config):
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


def test_split_and_save_data_ndarray(tmpdir, default_config):
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


def test_split_and_save_data_sparse(default_config, tmpdir):
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


def test_is_sanitizable():

    flag = is_sanitizable(("O"))
    assert flag

    flag = is_sanitizable(("c1ccccc1(C)(C)"))
    assert not flag


def test_reaction_to_fingerprint(default_config):
    default_config["fingerprint_len"] = 10
    product_smiles = "[Cl:1][c:2]1[c:3]([C:4](=[O:5])[C:12]([F:11])([F:13])[F:14])[cH:6][c:7]([F:10])[cH:8][cH:9]1"  # noqa
    reactants_smiles = "C[O:5][C:4](=O)[c:3]1[c:2]([Cl:1])[cH:9][cH:8][c:7]([F:10])[cH:6]1.C[Si](C)(C)[C:12]([F:11])([F:13])[F:14].COCCOC.Cl.[Cs+].[F-].[Na+].[OH-]"  # noqa

    fingerprint = reaction_to_fingerprints(
        [product_smiles, reactants_smiles], default_config
    )

    assert list(fingerprint) == [-1, -1, -1, 0, -1, 0, 0, -1, 0, -1]


def test_expansion_model_sequence_loading(expansion_model_sequence):

    assert expansion_model_sequence.input_dim == 10
    assert expansion_model_sequence.output_dim == 100


def test_expansion_model_sequence_slicing(expansion_model_sequence, default_config):
    seq = expansion_model_sequence

    xbatch, ybatch = seq[1]

    assert xbatch.shape[0] == default_config["batch_size"]
    assert ybatch.shape[0] == default_config["batch_size"]
    assert xbatch.shape[1] == expansion_model_sequence.input_dim
    assert ybatch.shape[1] == expansion_model_sequence.output_dim


def test_filter_model_sequence_loading(filter_model_sequence):

    assert filter_model_sequence.input_dim == 10


def test_filter_model_sequence_slicing(filter_model_sequence, default_config):
    seq = filter_model_sequence

    xbatch, ybatch = seq[1]

    assert xbatch[0].shape[0] == default_config["batch_size"]
    assert xbatch[1].shape[0] == default_config["batch_size"]
    assert ybatch.shape[0] == default_config["batch_size"]
    assert xbatch[0].shape[1] == filter_model_sequence.input_dim
    assert xbatch[1].shape[1] == filter_model_sequence.input_dim


def test_reactants_molecules():
    reactants_str = "C[O:5][C:4](=O)[c:3]1[c:2]([Cl:1])[cH:9][cH:8][c:7]([F:10])[cH:6]1.C[Si](C)(C)[C:12]([F:11])([F:13])[F:14].COCCOC.Cl.[Cs+].[F-]"  # noqa

    mols = create_reactants_molecules(reactants_str)

    assert len(mols) == 2
    expected_smiles = (
        "C[O:5][C:4](=O)[c:3]1[c:2]([Cl:1])[cH:9][cH:8][c:7]([F:10])[cH:6]1"
    )
    assert mols[0].smiles == expected_smiles
    assert mols[1].smiles == "C[Si](C)(C)[C:12]([F:11])([F:13])[F:14]"


def test_reverse_template():
    retro_template = "([C:2]-[C:3](=[O;D1;H0:4])-[N;H0;D3;+0:5](-[CH3;D1;+0:1])-[c:6])>>(I-[CH3;D1;+0:1]).([C:2]-[C:3](=[O;D1;H0:4])-[NH;D2;+0:5]-[c:6])"  # noqa
    expected = "(I-[CH3;D1;+0:1]).([C:2]-[C:3](=[O;D1;H0:4])-[NH;D2;+0:5]-[c:6])>>([C:2]-[C:3](=[O;D1;H0:4])-[N;H0;D3;+0:5](-[CH3;D1;+0:1])-[c:6])"  # noqa

    assert reverse_template(retro_template) == expected


def test_reaction_hash():
    reactants_str = "C[O:5][C:4](=O)[c:3]1[c:2]([Cl:1])[cH:9][cH:8][c:7]([F:10])[cH:6]1.C[Si](C)(C)[C:12]([F:11])([F:13])[F:14].COCCOC.Cl.[Cs+].[F-].[Na+].[OH-]"  # noqa
    product = Molecule(
        smiles="[Cl:1][c:2]1[c:3]([C:4](=[O:5])[C:12]([F:11])([F:13])[F:14])[cH:6][c:7]([F:10])[cH:8][cH:9]1"
    )
    expected = "c92d362cc28df004939d1af699cc8038a2f82ddcdf37eacf3c7366db"

    assert reaction_hash(reactants_str, product) == expected


@pytest.fixture
def library_data(shared_datadir, default_config):
    default_config["library_headers"] = [
        "index",
        "reaction_hash",
        "reactants",
        "products",
        "retro_template",
        "template_hash",
        "template_code",
    ]
    default_config["negative_data"]["random_trials"] = 10
    filename = str(shared_datadir / "make_false_template_library.csv")
    library = pd.read_csv(
        filename,
        index_col=False,
        header=None,
        names=default_config["library_headers"],
    )
    return library, default_config


def test_strict_application(library_data):
    library, config = library_data
    errors = []

    gen = strict_application(library, config, errors)
    new_df = next(gen)

    assert len(new_df) == 2
    assert not errors
    hashes = sorted(list(new_df.reaction_hash.values))
    assert hashes[0] == "a885c1d5fb93760b042b03a81d400d316908ce892f83a2de49e06042"
    assert hashes[1] == "fe2192395a637b200d94c29fddc4db38033456eacca0b8112052239c"

    assert next(gen) is None
    assert next(gen) is None

    with pytest.raises(StopIteration):
        next(gen)


def test_random_application(library_data):
    library, config = library_data
    errors = []

    gen = random_application(library, config, errors)

    assert next(gen) is None

    new_df = next(gen)

    assert len(new_df) == 1
    assert not errors
    hash_ = new_df.reaction_hash.values[0]
    assert hash_ == "acefe9ee01e5e4b6c84d4ed550f4819d791051f095d66abd60a80b25"

    assert next(gen) is None

    with pytest.raises(StopIteration):
        next(gen)


def test_recommender_application(library_data, mocker):
    mocked_load_model = mocker.patch(
        "aizynthfinder.training.make_false_products.load_keras_model"
    )
    mocked_model = mocked_load_model.return_value
    mocked_model.predict.return_value = np.asarray([0.5, 0.3, 0.2])
    library, config = library_data
    errors = []

    gen = recommender_application(library, config, errors)

    assert next(gen) is None

    new_df = next(gen)

    assert len(new_df) == 1
    assert not errors
    hash_ = new_df.reaction_hash.values[0]
    assert hash_ == "acefe9ee01e5e4b6c84d4ed550f4819d791051f095d66abd60a80b25"

    assert next(gen) is None

    with pytest.raises(StopIteration):
        next(gen)
