import os

import pytest
from PIL import Image, ImageDraw

from aizynthfinder.utils import image
from aizynthfinder.chem import TreeMolecule, Reaction


@pytest.fixture
def new_image():
    img = Image.new(mode="RGB", size=(300, 300), color="white")
    draw = ImageDraw.Draw(img)
    draw.point([5, 150, 100, 250], fill="black")
    return img


def test_crop(new_image):
    cropped = image.crop_image(new_image)

    assert cropped.width == 136
    assert cropped.height == 141
    assert cropped.getpixel((21, 21)) == (0, 0, 0)
    assert cropped.getpixel((116, 121)) == (0, 0, 0)


def test_rounded_rectangle(new_image):
    color = (255, 0, 0)
    modified = image.draw_rounded_rectangle(new_image, color=color)

    assert modified.getpixel((0, 150)) == color
    assert modified.getpixel((150, 0)) == color
    assert modified.getpixel((299, 150)) == color
    assert modified.getpixel((150, 299)) == color


@pytest.fixture
def setup_graphviz_graph():
    mol1 = TreeMolecule(smiles="CCCO", parent=None)
    reaction = Reaction(mol=mol1, smarts="")
    graph = image.GraphvizReactionGraph()

    graph.add_molecule(mol1, "green")
    graph.add_reaction(reaction)
    graph.add_edge(mol1, reaction)
    return graph


def test_graphviz_usage(mocker, tmpdir, setup_graphviz_graph):
    mkstemp_patch = mocker.patch("aizynthfinder.utils.image.tempfile.mkstemp")
    files = [
        (None, str(tmpdir / "graph1.dot")),
        (None, str(tmpdir / "graph2.dot")),
        (None, str(tmpdir / "img2.png")),
    ]
    mkstemp_patch.side_effect = files

    img = setup_graphviz_graph.to_image()

    assert img.height > 0
    assert img.width > 0
    for _, filename in files:
        assert os.path.exists(filename)


def test_graphviz_usage_exception_dot(mocker, tmpdir, setup_graphviz_graph):
    exists_patch = mocker.patch("aizynthfinder.utils.image.os.path.exists")
    exists_patch.return_value = False

    with pytest.raises(FileNotFoundError, match=".*'dot'.*"):
        setup_graphviz_graph.to_image()


def test_graphviz_usage_exception_neato(mocker, tmpdir, setup_graphviz_graph):
    exists_patch = mocker.patch("aizynthfinder.utils.image.os.path.exists")
    exists_patch.side_effect = [True, False]

    with pytest.raises(FileNotFoundError, match=".*'neato'.*"):
        setup_graphviz_graph.to_image()
