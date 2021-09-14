"""
Unit and regression test for the mmic_openff package.
"""

# Import package, test suite, and other packages as needed
import mmic_openff
import pytest
import sys
import os

import mmelemental as mm
from mmic_translator.models import schema_input_default
import mm_data
from .util import get_files

from openff.toolkit.topology import Molecule
from openff.toolkit.utils import get_data_file_path

molecules = [mfile for mfile in mm_data.mols if mfile.endswith(".sdf")]
molecules.extend(get_files("sdf", "molecules"))
molecules = [
    mfile
    for mfile in molecules
    if all(
        s not in mfile for s in ("multi", "conformer", "MiniDrug")
    )  # Can't handle multi-conformers for now ~ hackish
]


def test_mmic_openff_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_openff" in sys.modules


@pytest.mark.parametrize("mfile", molecules)
def test_mmic_to_mmschema(mfile: str, **kwargs):

    inputs = {
        "data_object": Molecule.from_file(mfile),
        "keywords": kwargs,
        "schema_version": 1,
        "schema_name": schema_input_default,
    }

    return mmic_openff.components.OpenFFToMolComponent.compute(inputs)


def test_mol_to_openff(**kwargs):
    mmol = mm.models.Molecule.from_file(get_data_file_path("molecules/toluene.pdb"))
    inputs = {
        "schema_object": mmol,
        "schema_version": 1,
        "schema_name": schema_input_default,
        "keywords": kwargs,
    }
    return mmic_openff.components.MolToOpenFFComponent.compute(inputs)


@pytest.mark.parametrize("mfile", molecules)
def test_io_methods(mfile: str, **kwargs):
    omol = mmic_openff.models.OpenFFMol.from_file(mfile)
    assert isinstance(omol.data, omol.dtype())

    omol.to_file("tmp.pdb")
    os.remove("tmp.pdb")

    mmol = omol.to_schema()
    assert isinstance(mmol, mm.models.Molecule)


@pytest.mark.parametrize("mfile", molecules)
def test_lossless_conv(mfile: str, **kwargs):
    pass
