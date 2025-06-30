""" This module contains a collection of routines to produce pretty images
"""

from __future__ import annotations

import atexit
import os
import shutil
import tempfile
from typing import TYPE_CHECKING

from jinja2 import Template
from rxnutils.routes.image import molecule_to_image, molecules_to_images

from aizynthfinder.chem import Molecule
from aizynthfinder.utils.paths import data_path

if TYPE_CHECKING:
    import networkx as nx

    from aizynthfinder.chem import FixedRetroReaction
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        PilColor,
        Sequence,
        Tuple,
        Union,
    )

IMAGE_FOLDER = tempfile.mkdtemp()


@atexit.register
def _clean_up_images() -> None:
    global IMAGE_FOLDER
    try:
        shutil.rmtree(IMAGE_FOLDER, ignore_errors=True)
    # pylint: disable=broad-except
    except Exception:  # noqa Don't care if we fail clean-up
        pass


def save_molecule_images(
    molecules: Sequence[Molecule], frame_colors: Sequence[PilColor], size: int = 300
) -> Dict[Molecule, str]:
    """
    Create images of a list of molecules and save them to disc
    a globally managed folder.

    :param molecules: the molecules to save as images
    :param frame_colors: the color of the frame around each image
    :param size: the sub-image size for each molecule
    :return: the filename of the created images
    """
    global IMAGE_FOLDER

    try:
        images = molecules_to_images(
            [mol.rd_mol for mol in molecules], frame_colors, size
        )
    # pylint: disable=broad-except
    except Exception:  # noqa
        images = [
            molecule_to_image(molecule.rd_mol, frame_color, size)
            for molecule, frame_color in zip(molecules, frame_colors)
        ]

    spec = {}
    for molecule, image_obj in zip(molecules, images):
        image_filepath = os.path.join(IMAGE_FOLDER, f"{molecule.inchi_key}.png")
        image_obj.save(image_filepath)
        spec[molecule] = image_filepath
    return spec


def make_visjs_page(
    filename: str,
    molecules: Sequence[Molecule],
    reactions: Sequence[FixedRetroReaction],
    edges: Union[Sequence[Tuple[Any, Any]], nx.digraph.OutEdgeView],
    frame_colors: Sequence[PilColor],
    hierarchical: bool = False,
) -> None:
    """
    Create HTML code of a bipartite graph of molecules and reactions
    using the vis.js network library.

    Package the created HTML page and all images as tar-ball.

    :param filename: the basename of the archive
    :param molecules: the molecules nodes
    :param reactions: the reaction nodes
    :param edges: the edges of the graph
    :param frame_colors: the color of the frame around each image
    :param hierarchical: if True, will produce a hierarchical layout
    """
    mol_spec = save_molecule_images(molecules, frame_colors)

    template_filepath = os.path.join(data_path(), "templates", "reaction_tree.thtml")
    with open(template_filepath, "r") as fileobj:
        template = Template(fileobj.read())
    template.globals["id"] = id  # type: ignore

    tmpdir = tempfile.mkdtemp()
    for image_filepath in mol_spec.values():
        shutil.copy(image_filepath, tmpdir)
    mol_spec = {molecule: os.path.basename(path) for molecule, path in mol_spec.items()}

    input_name = os.path.join(tmpdir, "route.html")
    with open(input_name, "w") as fileobj:
        fileobj.write(
            template.render(
                molecules=mol_spec,
                reactions=reactions,
                edges=edges,
                hierarchical=hierarchical,
            )
        )

    basename, _ = os.path.splitext(filename)
    shutil.make_archive(basename, "tar", root_dir=tmpdir)
