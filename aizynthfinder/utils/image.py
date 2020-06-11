""" This module contains a collection of routines to produce pretty images
"""
import sys
import subprocess
import os
import tempfile

from PIL import Image, ImageDraw
from rdkit.Chem import Draw
from rdkit import Chem


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


class GraphvizReactionGraph:
    """
    Class to create a reaction graph image using the graphviz set of tools

    The .dot files that are produced by this class is specifally made for this
    application. It is not at all intended for generic usage.
    """

    def __init__(self):
        self._graph_props = {"layout": "dot", "rankdir": "RL", "splines": "ortho"}
        self._mol_props = {"label": "", "color": "white", "shape": "none"}
        self._react_props = {
            "label": "",
            "fillcolor": "black",
            "shape": "circle",
            "style": "filled",
            "width": 0.1,
            "fixedsize": "true",
        }
        self._lines = [
            'strict digraph "" {',
            f"\t graph [{self._get_props_strings(self._graph_props)}\n\t];",
            'node [label="\\N"];',
        ]

    def add_edge(self, node1, node2):
        """ Add an edge to the graph between to nodes
        """
        self._lines.append(f'\t{id(node1)} -> {id(node2)} [arrowhead="none"];')

    def add_molecule(self, molecule, frame_color):
        """
        Add a node that is a molecule to the graph.
        The image of the molecule will have a rectangular frame around it with a given color.
        """
        img = molecule_to_image(molecule, frame_color=frame_color)
        _, filename = tempfile.mkstemp(suffix=".png")
        img.save(filename)
        self._mol_props["image"] = filename
        self._lines.append(
            f"\t{id(molecule)} [{self._get_props_strings(self._mol_props)}\n\t];"
        )

    def add_reaction(self, reaction):
        """ Add a node that is a reaction to the graph.
        """
        self._lines.append(
            f"\t{id(reaction)} [{self._get_props_strings(self._react_props)}\n\t];"
        )

    def to_image(self):
        """
        Produce the image using the 'dot' program for layout and
        'neato' for creating the image.

        :raises FileNotFoundError: if the intermediate files could not be created
        :return: the image of the reaction
        :rtype: PIL.Image
        """
        ext = ".bat" if sys.platform.startswith("win") else ""

        self._lines.append("}")
        _, input_name = tempfile.mkstemp(suffix=".dot")
        with open(input_name, "w") as fileobj:
            fileobj.write("\n".join(self._lines))

        _, output_dot = tempfile.mkstemp(suffix=".dot")
        subprocess.call([f"dot{ext}", "-T", "dot", f"-o{output_dot}", input_name])
        if not os.path.exists(output_dot):
            raise FileNotFoundError(
                "Could not produce graph with layout - check that 'dot' command is in path"
            )

        _, output_img = tempfile.mkstemp(suffix=".img")
        subprocess.call([f"neato{ext}", "-T", "png", f"-o{output_img}", output_dot])
        if not os.path.exists(output_img):
            raise FileNotFoundError(
                "Could not produce reaction image - check that 'neato' command is in path"
            )

        return Image.open(output_img)

    @staticmethod
    def _get_props_strings(props):
        return ",\n\t\t".join(f'{key}="{value}"' for key, value in props.items())
