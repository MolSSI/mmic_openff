"""
Unit and regression test for the mmic_openff package.
"""

# Import package, test suite, and other packages as needed
import mmic_openff
import pytest
import sys
import os

import mmelemental as mm
import mm_data

from openff.toolkit.topology import Molecule
from openff.toolkit.utils import get_data_file_path


def test_mmic_openff_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_openff" in sys.modules


def test_mmic_to_mol_from_sdf(**kwargs):
    sdf_file_path = get_data_file_path("molecules/ethanol.sdf")

    inputs = {
        "data_object": Molecule.from_file(sdf_file_path),
        "keywords": kwargs,
        "schema_version": 1,
        "schema_name": "mmschema",
    }

    return mmic_openff.components.OpenFFToMolComponent.compute(inputs)


def test_mol_to_openff(**kwargs):
    mmol = mm.models.Molecule.from_file(get_data_file_path("molecules/toluene.pdb"))
    inputs = {
        "schema_object": mmol,
        "schema_version": 1,
        "schema_name": "mmschema",
        "keywords": kwargs,
    }
    return mmic_openff.components.MolToOpenFFComponent.compute(inputs)


def test_io_methods(**kwargs):
    omol = mmic_openff.models.OpenFFMol.from_file(
        get_data_file_path("molecules/ethanol.sdf")
    )
    assert isinstance(omol.data, omol.dtype())

    omol.to_file("tmp.pdb")
    os.remove("tmp.pdb")
    # mmol = omol.to_schema()
    # assert isinstance(mmol, mm.models.molecule.Molecule)
