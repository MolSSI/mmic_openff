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

amber99sb = mm_data.ffs["amber99sb.xml"]


@pytest.mark.skip
def test_mmic_to_ff_from_xml(**kwargs):
    inputs = {
        "data_object": ForceField(amber99sb),
        "keywords": kwargs,
        "schema_version": 1,
        "schema_name": "my_schema",
    }
    return mmic_openff.components.OpenMMToFFComponent.compute(inputs)


@pytest.mark.skip("Cannot test this right now.")
def test_ff_to_openmm(**kwargs):
    ff = mm.models.ForceField.from_file(amber99sb)
    inputs = {"schema_object": ff, "keywords": kwargs}
    return mmic_openff.components.MolToOpenMMComponent.compute(inputs)


@pytest.mark.skip
def test_io_methods(**kwargs):
    ff = mmic_openff.models.OpenMMFF.from_file(amber99sb)
    assert isinstance(ff.data, ff.dtype())
