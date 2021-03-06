""" Module containing classes and routines used in training of policies.
"""
import os

import numpy as np
from tensorflow.keras.layers import Dense, Dropout, Input, Dot
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import Sequence
from tensorflow.keras.callbacks import (
    EarlyStopping,
    CSVLogger,
    ModelCheckpoint,
    ReduceLROnPlateau,
)
from tensorflow.keras import regularizers
from sklearn.utils import shuffle
from scipy import sparse

from aizynthfinder.utils.models import top10_acc, top50_acc


class _InMemorySequence(Sequence):
    def __init__(self, config, dataset_label):
        self.batch_size = config["batch_size"]
        input_filename = config.filename(dataset_label + "_inputs")
        label_filename = config.filename(dataset_label + "_labels")
        self.input_matrix = self._load_data(input_filename)
        self.label_matrix = self._load_data(label_filename)
        self.input_dim = self.input_matrix.shape[1]

    def __len__(self):
        return int(np.ceil(self.label_matrix.shape[0] / float(self.batch_size)))

    def _load_data(self, filename):
        try:
            return sparse.load_npz(filename)
        except ValueError:
            return np.load(filename)["arr_0"]

    def _make_slice(self, idx):
        if idx < 0 or idx >= len(self):
            raise IndexError("index out of range")

        start = idx * self.batch_size
        end = (idx + 1) * self.batch_size
        return slice(start, end)


class ExpansionModelSequence(_InMemorySequence):
    """
    Custom sequence class to keep sparse, pre-computed matrices in memory.
    Batches are created dynamically by slicing the in-memory arrays
    The data will be shuffled on each epoch end

    :ivar batch_size: the batch size as read from config
    :vartype batch_size: int
    :ivar input_matrix: the loaded input matrix
    :vartype input_matrix: scipy.sparse object
    :ivar label_matrix: the loaded label matrix
    :vartype label_matrix: scipy.sparse object
    :ivar input_dim: the input size (fingerprint size)
    :vartype input_dim: int
    :ivar output_dim: the output size (number of templates)
    :vartype output_dim: int

    :param config: the settings
    :type config: Config
    :param dataset_label: the label of set, e.g. training, testing or validation
    :type dataset_label: str
    """

    def __init__(self, config, dataset_label):
        super().__init__(config, dataset_label)
        self.output_dim = self.label_matrix.shape[1]

    def __getitem__(self, idx):
        idx = self._make_slice(idx)
        return self.input_matrix[idx].toarray(), self.label_matrix[idx].toarray()

    def on_epoch_end(self):
        self.input_matrix, self.label_matrix = shuffle(
            self.input_matrix, self.label_matrix, random_state=0
        )


class FilterModelSequence(_InMemorySequence):
    """
    Custom sequence class to keep sparse, pre-computed matrices in memory.
    Batches are created dynamically by slicing the in-memory arrays
    The data will be shuffled on each epoch end

    :ivar batch_size: the batch size as read from config
    :vartype batch_size: int
    :ivar input_matrix: the loaded input matrix for the products
    :vartype input_matrix: scipy.sparse object
    :ivar input_matrix2: the loaded input matrix for the reactions
    :vartype input_matrix2: scipy.sparse object
    :ivar label_matrix: the loaded label matrix
    :vartype label_matrix: np.ndarray
    :ivar input_dim: the input size (fingerprint size)
    :vartype input_dim: int

    :param config: the settings
    :type config: Config
    :param dataset_label: the label of set, e.g. training, testing or validation
    :type dataset_label: str
    """

    def __init__(self, config, dataset_label):
        super().__init__(config, dataset_label)
        filename = config.filename(dataset_label + "_inputs2")
        self.input_matrix2 = self._load_data(filename)

    def __getitem__(self, idx):
        idx = self._make_slice(idx)
        return (
            [self.input_matrix[idx].toarray(), self.input_matrix2[idx].toarray()],
            self.label_matrix[idx],
        )

    def on_epoch_end(self):
        self.input_matrix, self.input_matrix2, self.label_matrix = shuffle(
            self.input_matrix, self.input_matrix2, self.label_matrix, random_state=0
        )


def _setup_callbacks(config):
    early_stopping = EarlyStopping(monitor="val_loss", patience=10)
    csv_logger = CSVLogger(config.filename("_keras_training.log"), append=True)

    checkpoint_path = os.path.join(config["output_path"], "checkpoints")
    if not os.path.exists(checkpoint_path):
        os.mkdir(checkpoint_path)
    checkpoint = ModelCheckpoint(
        os.path.join(checkpoint_path, "keras_model.hdf5"),
        monitor="loss",
        save_best_only=True,
    )

    reduce_lr = ReduceLROnPlateau(
        monitor="val_loss",
        factor=0.5,
        patience=5,
        verbose=0,
        mode="auto",
        min_delta=0.000001,
        cooldown=0,
        min_lr=0,
    )
    return [early_stopping, csv_logger, checkpoint, reduce_lr]


def _train_keras_model(model, train_seq, valid_seq, loss, metrics, config):
    adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0)

    model.compile(
        optimizer=adam, loss=loss, metrics=metrics,
    )

    model.fit_generator(
        train_seq,
        steps_per_epoch=None,
        epochs=config["epochs"],
        verbose=1,
        callbacks=_setup_callbacks(config),
        validation_data=valid_seq,
        validation_steps=None,
        class_weight=None,
        max_queue_size=20,
        workers=20,
        use_multiprocessing=False,
        shuffle=True,
        initial_epoch=0,
    )


def train_expansion_keras_model(config):
    """
    Train a expansion policy

    :param config: the settings
    :type config: Config
    """
    train_seq = ExpansionModelSequence(config, "training")
    valid_seq = ExpansionModelSequence(config, "validation")

    model = Sequential()
    model.add(
        Dense(
            config["model"]["hidden_nodes"],
            input_shape=(train_seq.input_dim,),
            activation="elu",
            kernel_regularizer=regularizers.l2(0.001),
        )
    )
    model.add(Dropout(config["model"]["drop_out"]))
    model.add(Dense(train_seq.output_dim, activation="softmax"))

    _train_keras_model(
        model,
        train_seq,
        valid_seq,
        "categorical_crossentropy",
        ["accuracy", "top_k_categorical_accuracy", top10_acc, top50_acc],
        config,
    )


def train_filter_keras_model(config):

    train_seq = FilterModelSequence(config, "training")
    valid_seq = FilterModelSequence(config, "validation")

    product_input_layer = Input(shape=(config["fingerprint_len"],))
    product_dense_layer = Dense(config["model"]["hidden_nodes"], activation="elu")(
        product_input_layer
    )
    product_droput_layer = Dropout(config["model"]["drop_out"])(product_dense_layer)
    reaction_input_layer = Input(shape=(config["fingerprint_len"],))
    reaction_dense_layer = Dense(config["model"]["hidden_nodes"], activation="elu")(
        reaction_input_layer
    )
    cosine_layer = Dot(-1, normalize=True)([product_droput_layer, reaction_dense_layer])
    output_layer = Dense(1, activation="sigmoid")(cosine_layer)
    model = Model(
        inputs=[product_input_layer, reaction_input_layer], outputs=output_layer
    )

    _train_keras_model(
        model, train_seq, valid_seq, "binary_crossentropy", ["accuracy"], config
    )


def train_recommender_keras_model(config):

    train_seq = ExpansionModelSequence(config, "training")
    valid_seq = ExpansionModelSequence(config, "validation")

    model = Sequential()
    model.add(
        Dense(
            config["model"]["hidden_nodes"],
            input_shape=(config["fingerprint_len"],),
            activation="elu",
        )
    )
    model.add(Dropout(config["model"]["drop_out"]))
    model.add(Dense(train_seq.output_dim, activation="softmax"))

    _train_keras_model(
        model,
        train_seq,
        valid_seq,
        "categorical_crossentropy",
        ["accuracy", "top_k_categorical_accuracy", top10_acc, top50_acc],
        config,
    )
