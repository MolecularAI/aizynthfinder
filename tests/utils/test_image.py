import os
from tarfile import TarFile

import pytest
from PIL import Image, ImageDraw

from aizynthfinder.utils import image
from aizynthfinder.chem import TreeMolecule, TemplatedRetroReaction


@pytest.fixture
def new_image():
    img = Image.new(mode="RGB", size=(300, 300), color="white")
    draw = ImageDraw.Draw(img)
    draw.point([5, 150, 100, 250], fill="black")
    return img


@pytest.fixture
def setup_graphviz_graph():
    mol1 = TreeMolecule(smiles="CCCO", parent=None)
    reaction = TemplatedRetroReaction(mol=mol1, smarts="")

    return [mol1], [reaction], [(mol1, reaction)], ["green"]


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


def test_save_molecule_images():
    nfiles = len(os.listdir(image.IMAGE_FOLDER))

    mols = [
        TreeMolecule(smiles="CCCO", parent=None),
        TreeMolecule(smiles="CCCO", parent=None),
        TreeMolecule(smiles="CCCCO", parent=None),
    ]

    image.save_molecule_images(mols, ["green", "green", "green"])

    assert len(os.listdir(image.IMAGE_FOLDER)) == nfiles + 2

    image.save_molecule_images(mols, ["green", "orange", "green"])

    assert len(os.listdir(image.IMAGE_FOLDER)) == nfiles + 2


def test_graphviz_usage(mocker, tmpdir, setup_graphviz_graph):
    mkstemp_patch = mocker.patch("aizynthfinder.utils.image.tempfile.mkstemp")
    files = [
        (None, str(tmpdir / "graph1.dot")),
        (None, str(tmpdir / "img2.png")),
    ]
    mkstemp_patch.side_effect = files
    molecules, reactions, edges, frame_colors = setup_graphviz_graph

    img = image.make_graphviz_image(molecules, reactions, edges, frame_colors)

    assert img.height > 0
    assert img.width > 0
    for _, filename in files:
        assert os.path.exists(filename)


def test_graphviz_usage_exception_dot(mocker, tmpdir, setup_graphviz_graph):
    exists_patch = mocker.patch("aizynthfinder.utils.image.os.path.exists")
    exists_patch.side_effect = [False, True]
    molecules, reactions, edges, frame_colors = setup_graphviz_graph

    img = image.make_graphviz_image(molecules, reactions, edges, frame_colors)
    assert img.height > 0
    assert img.width > 0


def test_graphviz_usage_exception_dot_both(mocker, tmpdir, setup_graphviz_graph):
    exists_patch = mocker.patch("aizynthfinder.utils.image.os.path.exists")
    exists_patch.return_value = False
    molecules, reactions, edges, frame_colors = setup_graphviz_graph

    with pytest.raises(FileNotFoundError, match=".*'dot'.*"):
        image.make_graphviz_image(molecules, reactions, edges, frame_colors)


def test_visjs_page(mocker, tmpdir, setup_graphviz_graph):
    mkdtemp_patch = mocker.patch("aizynthfinder.utils.image.tempfile.mkdtemp")
    mkdtemp_patch.return_value = str(tmpdir / "tmp")
    os.mkdir(tmpdir / "tmp")
    molecules, reactions, edges, frame_colors = setup_graphviz_graph
    filename = str(tmpdir / "arch.tar")

    image.make_visjs_page(filename, molecules, reactions, edges, frame_colors)

    assert os.path.exists(filename)
    with TarFile(filename) as tarobj:
        assert "./route.html" in tarobj.getnames()
        assert len([name for name in tarobj.getnames() if name.endswith(".png")]) == 1
