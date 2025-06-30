import os
from tarfile import TarFile
import pytest

from aizynthfinder.utils import image
from aizynthfinder.chem import TreeMolecule, TemplatedRetroReaction


@pytest.fixture
def setup_graph():
    mol1 = TreeMolecule(smiles="CCCO", parent=None)
    reaction = TemplatedRetroReaction(mol=mol1, smarts="")

    return [mol1], [reaction], [(mol1, reaction)], ["green"]


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
