import os
import shutil
import json
from tarfile import TarFile
from pathlib import Path

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
def setup_graph():
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


def test_visjs_page(mocker, tmpdir, setup_graph):
    mkdtemp_patch = mocker.patch("aizynthfinder.utils.image.tempfile.mkdtemp")
    mkdtemp_patch.return_value = str(tmpdir / "tmp")
    os.mkdir(tmpdir / "tmp")
    molecules, reactions, edges, frame_colors = setup_graph
    filename = str(tmpdir / "arch.tar")

    image.make_visjs_page(filename, molecules, reactions, edges, frame_colors)

    assert os.path.exists(filename)
    with TarFile(filename) as tarobj:
        assert "./route.html" in tarobj.getnames()
        assert len([name for name in tarobj.getnames() if name.endswith(".png")]) == 1


def test_image_factory(request):
    route_path = Path(request.fspath).parent.parent / "data" / "branched_route.json"
    with open(route_path, "r") as fileobj:
        dict_ = json.load(fileobj)
    dict_["children"][0]["children"][1]["hide"] = True

    factory0 = image.RouteImageFactory(dict_)

    factory_tighter = image.RouteImageFactory(dict_, margin=50)
    assert factory0.image.width == factory_tighter.image.width + 150
    assert factory0.image.height == factory_tighter.image.height + 175

    factory_hidden = image.RouteImageFactory(dict_, show_all=False)
    assert factory0.image.width == factory_hidden.image.width
    assert factory0.image.height > factory_hidden.image.height
