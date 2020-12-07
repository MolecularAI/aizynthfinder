""" This module contains a collection of routines to produce pretty images
"""
import sys
import subprocess
import os
import tempfile
import atexit
import shutil

from jinja2 import Template
from PIL import Image, ImageDraw
from rdkit.Chem import Draw
from rdkit import Chem

from aizynthfinder.utils.paths import data_path

IMAGE_FOLDER = tempfile.mkdtemp()


@atexit.register
def _clean_up_images():
    global IMAGE_FOLDER
    try:
        shutil.rmtree(IMAGE_FOLDER, ignore_errors=True)
    except Exception:  # Don't care if we fail clean-up
        pass


def molecule_to_image(mol, frame_color):
    """
    Create a pretty image of a molecule,
    with a colored frame around it

    :param mol: the molecule
    :type mol: Molecule
    :param frame_color: the color of the frame
    :type frame_color: tuple of int or str
    :return: the produced image
    :rtype: PIL image
    """
    mol = Chem.MolFromSmiles(mol.smiles)
    img = Draw.MolToImage(mol)
    cropped_img = crop_image(img)
    return draw_rounded_rectangle(cropped_img, frame_color)


def crop_image(img, margin=20):
    """
    Crop an image by removing white space around it

    :param img: the image to crop
    :type img: PIL image
    :param margin: padding, defaults to 20
    :type margin: int, optional
    :return: the cropped image
    :rtype: PIL image
    """
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


def draw_rounded_rectangle(img, color, arc_size=20):
    """
    Draw a rounded rectangle around an image

    :param img: the image to draw upon
    :type img: PIL image
    :param color: the color of the rectangle
    :type color: tuple or str
    :param arc_size: the size of the corner, defaults to 20
    :type arc_size: int, optional
    :return: the new image
    :rtype: PIL image
    """
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


def save_molecule_images(molecules, frame_colors):
    """
    Create images of a list of molecules and save them to disc
    a globally managed folder.

    :param molecules: the molecules to save as images
    :type molecules: list of Molecule
    :param frame_colors: the color of the frame around each image
    :type frame_colors: list of str
    :return: the filename of the created images
    :rtype: dict
    """
    global IMAGE_FOLDER
    spec = {}
    for molecule, frame_color in zip(molecules, frame_colors):
        image_filepath = os.path.join(
            IMAGE_FOLDER, f"{molecule.inchi_key}_{frame_color}.png"
        )
        if not os.path.exists(image_filepath):
            image_obj = molecule_to_image(molecule, frame_color)
            image_obj.save(image_filepath)
        spec[molecule] = image_filepath
    return spec


def make_graphviz_image(molecules, reactions, edges, frame_colors):
    """
    Create an image of a bipartite graph of molecules and reactions
    using the dot program of graphviz

    :param molecules: the molecules nodes
    :type molecules: list of Molecules
    :param reactions: the reaction nodes
    :type reactions: list of Reactions
    :param edges: the edges of the graph
    :type edges: list of tuples
    :param frame_colors: the color of the frame around each image
    :type frame_colors: list of str
    :raises FileNotFoundError: if the image could not be produced
    :return: the create image
    :rtype: PIL.Image
    """

    def _create_image(use_splines):
        txt = template.render(
            molecules=mol_spec,
            reactions=reactions,
            edges=edges,
            use_splines=use_splines,
        )
        _, input_name = tempfile.mkstemp(suffix=".dot")
        with open(input_name, "w") as fileobj:
            fileobj.write(txt)

        _, output_img = tempfile.mkstemp(suffix=".png")
        ext = ".bat" if sys.platform.startswith("win") else ""
        subprocess.call([f"dot{ext}", "-T", "png", f"-o{output_img}", input_name])
        if not os.path.exists(output_img) or os.path.getsize(output_img) == 0:
            raise FileNotFoundError(
                "Could not produce graph with layout - check that 'dot' command is in path"
            )
        return output_img

    mol_spec = save_molecule_images(molecules, frame_colors)

    template_filepath = os.path.join(data_path(), "templates", "reaction_tree.dot")
    with open(template_filepath, "r") as fileobj:
        template = Template(fileobj.read())
    template.globals["id"] = id

    try:
        output_img = _create_image(use_splines=True)
    except FileNotFoundError:
        output_img = _create_image(use_splines=False)

    return Image.open(output_img)


def make_visjs_page(
    filename, molecules, reactions, edges, frame_colors, hierarchical=False
):
    """
    Create HTML code of a bipartite graph of molecules and reactions
    using the vis.js network library.

    Package the created HTML page and all images as tar-ball.

    :param filename: the basename of the archive
    :type filename: str
    :param molecules: the molecules nodes
    :type molecules: list of Molecules
    :param reactions: the reaction nodes
    :type reactions: list of Reactions
    :param edges: the edges of the graph
    :type edges: list of tuples
    :param frame_colors: the color of the frame around each image
    :type frame_colors: list of str
    :param hierarchical: if True, will produce a hierarchical layout
    :type hierarchical: bool, optional
    """
    mol_spec = save_molecule_images(molecules, frame_colors)

    template_filepath = os.path.join(data_path(), "templates", "reaction_tree.thtml")
    with open(template_filepath, "r") as fileobj:
        template = Template(fileobj.read())
    template.globals["id"] = id

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
