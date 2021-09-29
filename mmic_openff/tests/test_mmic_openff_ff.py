"""
Unit and regression test for the mmic_openff package.
"""

# Import package, test suite, and other packages as needed
import mmic_openff
from openff.toolkit.typing.engines.smirnoff import ForceField
import pytest
import sys
import os
import mmelemental as mm
import mm_data

forcefields = [
    "openff-1.2.0.offxml",
    "openff-2.0.0.offxml",
    pytest.param(mm_data.ffs["amber99sb.xml"], marks=pytest.mark.skip),
]


@pytest.mark.parametrize("ff", forcefields)
def test_mmic_openff_from_xml(ff, **kwargs):
    inputs = {
        "data_object": ForceField(ff),
        "keywords": kwargs,
        "schema_version": mmic_openff.mmic_openff._mmschema_max_version,
        "schema_name": mm.models.ForceField.default_schema_name,
    }
    return mmic_openff.components.OpenFFToFFComponent.compute(inputs)


@pytest.mark.parametrize("ff", forcefields)
def test_mmic_openff_from_mmschema(ff, **kwargs):
    mmff = test_mmic_openff_from_xml(ff).schema_object
    mmff = mm.models.ForceField(**mmff.dict())
    inputs = {
        "schema_object": mmff,
        "schema_name": mmff.schema_name,
        "schema_version": mmff.schema_version,
        "keywords": kwargs,
    }
    return mmic_openff.components.FFToOpenFFComponent.compute(inputs)


@pytest.mark.skip
def test_io_methods(**kwargs):
    ff = mmic_openff.models.OpenMMFF.from_file(amber99sb)
    assert isinstance(ff.data, ff.dtype)
