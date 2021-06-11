""" Module routines for creating negative data for the training for filter policies
"""
from __future__ import annotations
import random
import argparse
from typing import TYPE_CHECKING

import pandas as pd
import tqdm

import aizynthfinder.utils.logging  # pylint: disable=unused-import
from aizynthfinder.chem import Molecule, Reaction, MoleculeException
from aizynthfinder.training.utils import (
    Config,
    create_reactants_molecules,
    reverse_template,
    reaction_hash,
    reactants_to_fingerprint,
)
from aizynthfinder.utils.models import CUSTOM_OBJECTS, load_keras_model

if TYPE_CHECKING:
    from aizynthfinder.utils.type_utils import (
        Iterable,
        Optional,
        Tuple,
        Any,
        Callable,
        List,
    )

    _DfGenerator = Iterable[Optional[pd.DataFrame]]


class _ReactionException(Exception):
    pass


def random_application(library: pd.DataFrame, config: Config, _) -> _DfGenerator:
    """
    Apply a random template to each row of reactants in the library
    to make new, false reactions.

    :param library: the reaction library
    :param config: the configuration
    :yield: a new DataFrame with a false reaction for each row if a match could be found, otherwise None
    """
    random.seed(100)  # To have reproducible results

    def random_sampler(_):
        for _ in range(config["negative_data"]["random_trials"]):
            template_index = random.randrange(len(library))
            yield library.iloc[template_index]

    yield from _sample_library(library, config, random_sampler)


def recommender_application(library: pd.DataFrame, config: Config, _) -> _DfGenerator:
    """
    Apply a template recommended by a neural network
    to each row of reactants in the library in order to make new, false reactions.

    :param library: the reaction library
    :param config: the configuration
    :yield: a new DataFrame with a false reaction for each row if a match could be found, otherwise None
    """

    model = load_keras_model(
        config["negative_data"]["recommender_model"], custom_objects=CUSTOM_OBJECTS
    )
    topn = config["negative_data"]["recommender_topn"]

    def prediction_sampler(row):
        fingerprint = reactants_to_fingerprint([row.reactants], config)
        fingerprint = fingerprint.reshape([1, config["fingerprint_len"]])
        prediction = model.predict(fingerprint).flatten()
        prediction_indices = prediction.argsort()[::-1][:topn]
        for template_code in prediction_indices:
            yield library[library.template_code == template_code].iloc[0]

    yield from _sample_library(library, config, prediction_sampler)


def strict_application(
    library: pd.DataFrame, config: Config, errors: list
) -> _DfGenerator:
    """
    Apply the recorded template to each row of reactants in the library
    to make new, false reactions.

    :param library: the reaction library
    :param config: the configuration
    :param errors: a list to fill with strings of produced errors
    :yield: a new DataFrame with a false reaction for each row if a match could be found, otherwise None
    """
    for _, row in library.iterrows():
        try:
            new_df = _apply_forward_reaction(row, config)
        except _ReactionException as err:
            errors.append(str(err))
        else:
            yield new_df


def _apply_forward_reaction(
    template_row: pd.Series, config: Config
) -> Optional[pd.DataFrame]:
    smarts_fwd = reverse_template(template_row.retro_template)
    mols = create_reactants_molecules(template_row.reactants)

    try:
        ref_mol = Molecule(smiles=template_row.products, sanitize=True)
    except MoleculeException as err:
        raise _ReactionException(
            f"reaction {template_row.reaction_hash} failed with msg {str(err)}"
        )

    try:
        products = Reaction(mols=mols, smarts=smarts_fwd).apply()
    except ValueError as err:
        raise _ReactionException(
            f"reaction {template_row.reaction_hash} failed with msg {str(err)}"
        )

    new_products = {product[0] for product in products if product[0] != ref_mol}
    if not new_products:
        return None

    correct_products = {
        product[0] for product in products if product[0].basic_compare(ref_mol)
    }
    if not correct_products:
        raise _ReactionException(
            f"reaction {template_row.reaction_hash} failed to produce correct product"
        )

    return _new_dataframe(
        template_row,
        config,
        nrows=len(new_products),
        reaction_hash=[
            reaction_hash(template_row.reactants, product) for product in new_products
        ],
        products=[product.smiles for product in new_products],
    )


def _get_config() -> Tuple[Config, str]:
    parser = argparse.ArgumentParser("Tool to generate artificial negative reactions")
    parser.add_argument("config", help="the filename to a configuration file")
    parser.add_argument(
        "method",
        choices=["strict", "random", "recommender"],
        help="the method to create random data",
    )
    args = parser.parse_args()

    return Config(args.config), args.method


def _new_dataframe(
    original: pd.Series, config: Config, nrows: int = 1, **kwargs: Any
) -> pd.DataFrame:
    dict_ = {"index": 0}
    for column in config["library_headers"][1:]:
        dict_[column] = kwargs.get(column, [original[column]] * nrows)
    return pd.DataFrame(dict_)


def _sample_library(
    library: pd.DataFrame, config: Config, sampler_func: Callable
) -> _DfGenerator:
    for _, row in library.iterrows():
        mols = create_reactants_molecules(row.reactants)
        try:
            ref_mol = Molecule(smiles=row.products, sanitize=True)
        except MoleculeException:
            yield None
            continue

        new_product = None
        for template_row in sampler_func(row):
            if row.template_hash == template_row.template_hash:
                continue
            smarts_fwd = reverse_template(template_row.retro_template)
            try:
                new_product = Reaction(mols=mols, smarts=smarts_fwd).apply()[0][0]
            except (ValueError, IndexError):
                continue
            if new_product.basic_compare(ref_mol):
                continue
            break  # If we have reached here, we have found a match that fits all criteria

        if not new_product:
            yield None
            continue

        # pylint: disable=undefined-loop-variable
        yield _new_dataframe(
            row,
            config,
            reaction_hash=[reaction_hash(row.reactants, new_product)],
            products=[new_product.smiles],
            classification=[""],
            retro_template=[template_row.retro_template],
            template_hash=[template_row.template_hash],
            selectivity=[0],
            outcomes=[1],
            template_code=[template_row.template_code],
        )


def main() -> None:
    """Entry-point for the make_false_products tool"""
    methods = {
        "strict": strict_application,
        "random": random_application,
        "recommender": recommender_application,
    }

    config, selected_method = _get_config()
    filename = config.filename("library")
    library = pd.read_csv(
        filename,
        index_col=False,
        header=None,
        names=config["library_headers"],
    )
    false_lib = pd.DataFrame({column: [] for column in config["library_headers"]})

    progress_bar = tqdm.tqdm(total=len(library))
    errors: List[str] = []
    for new_df in methods[selected_method](library, config, errors):
        if new_df is not None:
            false_lib = false_lib.append(new_df, sort=False)
        progress_bar.update(1)
    progress_bar.close()

    false_lib.to_csv(
        config.filename("false_library"),
        mode="w",
        header=False,
        index=False,
    )
    with open(config.filename("_errors.txt"), "w") as fileobj:
        fileobj.write("\n".join(errors))


if __name__ == "__main__":
    main()
