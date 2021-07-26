""" This module contains a collection of routines to produce pretty images
"""
from __future__ import annotations
import sys
import io
import subprocess
import os
import tempfile
import atexit
import shutil
from typing import TYPE_CHECKING

from jinja2 import Template
from PIL import Image, ImageDraw
from rdkit.Chem import Draw
from rdkit import Chem

from aizynthfinder.utils.paths import data_path

if TYPE_CHECKING:
    import networkx as nx

    # pylint: disable=ungrouped-imports
    from aizynthfinder.utils.type_utils import (
        Tuple,
        Any,
        Dict,
        Union,
        Sequence,
        PilImage,
        PilColor,
        List,
    )
    from aizynthfinder.chem import (
        Molecule,
        UniqueMolecule,
        RetroReaction,
        FixedRetroReaction,
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


def molecule_to_image(mol: Molecule, frame_color: PilColor) -> PilImage:
    """
    Create a pretty image of a molecule,
    with a colored frame around it

    :param mol: the molecule
    :param frame_color: the color of the frame
    :return: the produced image
    """
    mol = Chem.MolFromSmiles(mol.smiles)
    img = Draw.MolToImage(mol)
    cropped_img = crop_image(img)
    return draw_rounded_rectangle(cropped_img, frame_color)


def molecules_to_images(
    mols: Sequence[Molecule], frame_colors: Sequence[PilColor]
) -> List[PilImage]:
    """
    Create pretty images of molecules with a colored frame around each one of them.

    The molecules will be resized to be of similar sizes.

    :param mols: the molecules
    :param frame_colors: the color of the frame for each molecule
    :return: the produced images
    """
    size = 300
    # Make sanitized copies of all molecules
    mol_copies = [mol.make_unique() for mol in mols]
    for mol in mol_copies:
        mol.sanitize()

    all_mols = Draw.MolsToGridImage(
        [mol.rd_mol for mol in mol_copies],
        molsPerRow=len(mols),
        subImgSize=(size, size),
    )
    if not hasattr(all_mols, "crop"):  # Is not a PIL image
        fileobj = io.BytesIO(all_mols.data)
        all_mols = Image.open(fileobj)

    images = []
    for idx, frame_color in enumerate(frame_colors):
        image_obj = all_mols.crop((size * idx, 0, size * (idx + 1), size))
        image_obj = crop_image(image_obj)
        images.append(draw_rounded_rectangle(image_obj, frame_color))
    return images


def crop_image(img: PilImage, margin: int = 20) -> PilImage:
    """
    Crop an image by removing white space around it

    :param img: the image to crop
    :param margin: padding, defaults to 20
    :return: the cropped image
    """
    # pylint: disable=invalid-name
    # First find the boundaries of the white area
    x0_lim = img.width
    y0_lim = img.height
    x1_lim = 0
    y1_lim = 0
    for x in range(0, img.width):
        for y in range(0, img.height):
            if img.getpixel((x, y)) != (255, 255, 255):
                if x < x0_lim:
                    x0_lim = x
                if x > x1_lim:
                    x1_lim = x
                if y < y0_lim:
                    y0_lim = y
                if y > y1_lim:
                    y1_lim = y
    x0_lim = max(x0_lim, 0)
    y0_lim = max(y0_lim, 0)
    x1_lim = min(x1_lim + 1, img.width)
    y1_lim = min(y1_lim + 1, img.height)
    # Then crop to this area
    cropped = img.crop((x0_lim, y0_lim, x1_lim, y1_lim))
    # Then create a new image with the desired padding
    out = Image.new(
        img.mode,
        (cropped.width + 2 * margin, cropped.height + 2 * margin),
        color="white",
    )
    out.paste(cropped, (margin + 1, margin + 1))
    return out


def draw_rounded_rectangle(
    img: PilImage, color: PilColor, arc_size: int = 20
) -> PilImage:
    """
    Draw a rounded rectangle around an image

    :param img: the image to draw upon
    :param color: the color of the rectangle
    :param arc_size: the size of the corner, defaults to 20
    :return: the new image
    """
    # pylint: disable=invalid-name
    x0, y0, x1, y1 = img.getbbox()
    x1 -= 1
    y1 -= 1
    copy = img.copy()
    draw = ImageDraw.Draw(copy)
    arc_size_half = arc_size // 2
    draw.arc((x0, y0, arc_size, arc_size), start=180, end=270, fill=color)
    draw.arc((x1 - arc_size, y0, x1, arc_size), start=270, end=0, fill=color)
    draw.arc((x1 - arc_size, y1 - arc_size, x1, y1), start=0, end=90, fill=color)
    draw.arc((x0, y1 - arc_size, arc_size, y1), start=90, end=180, fill=color)
    draw.line((x0 + arc_size_half, y0, x1 - arc_size_half, y0), fill=color)
    draw.line((x1, arc_size_half, x1, y1 - arc_size_half), fill=color)
    draw.line((arc_size_half, y1, x1 - arc_size_half, y1), fill=color)
    draw.line((x0, arc_size_half, x0, y1 - arc_size_half), fill=color)
    return copy


def save_molecule_images(
    molecules: Sequence[Molecule], frame_colors: Sequence[PilColor]
) -> Dict[Molecule, str]:
    """
    Create images of a list of molecules and save them to disc
    a globally managed folder.

    :param molecules: the molecules to save as images
    :param frame_colors: the color of the frame around each image
    :return: the filename of the created images
    """
    global IMAGE_FOLDER

    try:
        images = molecules_to_images(molecules, frame_colors)
    except Exception:  # noqa
        images = [
            molecule_to_image(molecule, frame_color)
            for molecule, frame_color in zip(molecules, frame_colors)
        ]

    spec = {}
    for molecule, image_obj in zip(molecules, images):
        image_filepath = os.path.join(IMAGE_FOLDER, f"{molecule.inchi_key}.png")
        image_obj.save(image_filepath)
        spec[molecule] = image_filepath
    return spec


def make_graphviz_image(
    molecules: Union[Sequence[Molecule], Sequence[UniqueMolecule]],
    reactions: Union[Sequence[RetroReaction], Sequence[FixedRetroReaction]],
    edges: Sequence[Tuple[Any, Any]],
    frame_colors: Sequence[PilColor],
) -> PilImage:
    """
    Create an image of a bipartite graph of molecules and reactions
    using the dot program of graphviz

    :param molecules: the molecules nodes
    :param reactions: the reaction nodes
    :param edges: the edges of the graph
    :param frame_colors: the color of the frame around each image
    :raises FileNotFoundError: if the image could not be produced
    :return: the create image
    """

    def _create_image(use_splines):
        txt = template.render(
            molecules=mol_spec,
            reactions=reactions,
            edges=edges,
            use_splines=use_splines,
        )
        _, input_name = tempfile.mkstemp(suffix=".dot")
        with open(input_name, "w") as this_fileobj:
            this_fileobj.write(txt)

        _, output_img2 = tempfile.mkstemp(suffix=".png")
        ext = ".bat" if sys.platform.startswith("win") else ""
        subprocess.call([f"dot{ext}", "-T", "png", f"-o{output_img2}", input_name])
        if not os.path.exists(output_img2) or os.path.getsize(output_img2) == 0:
            raise FileNotFoundError(
                "Could not produce graph with layout - check that 'dot' command is in path"
            )
        return output_img2

    mol_spec = save_molecule_images(molecules, frame_colors)

    template_filepath = os.path.join(data_path(), "templates", "reaction_tree.dot")
    with open(template_filepath, "r") as fileobj:
        template = Template(fileobj.read())
    template.globals["id"] = id  # type: ignore

    try:
        output_img = _create_image(use_splines=True)
    except FileNotFoundError:
        output_img = _create_image(use_splines=False)

    return Image.open(output_img)


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
