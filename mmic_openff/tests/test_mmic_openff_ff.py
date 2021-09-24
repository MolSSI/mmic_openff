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
    "openff-2.0.0.offxml",
    #    pytest.param(mm_data.ffs["amber99sb.xml"], marks=pytest.mark.skip),
]


#@pytest.mark.parametrize("ff", forcefields)
@pytest.mark.skip("Skip temporarily.")
def test_mmic_to_ff_from_xml(ff, **kwargs):
    inputs = {
        "data_object": ForceField(ff),
        "keywords": kwargs,
        "schema_version": 1,
        "schema_name": mm.models.ForceField.default_schema_name,
    }
    return mmic_openff.components.OpenFFToFFComponent.compute(inputs)


@pytest.mark.skip("Cannot test this right now.")
def test_ff_to_openmm(**kwargs):
    ff = mm.models.ForceField.from_file(amber99sb)
    inputs = {"schema_object": ff, "keywords": kwargs}
    return mmic_openff.components.FFToOpenFFComponent.compute(inputs)


@pytest.mark.skip
def test_io_methods(**kwargs):
    ff = mmic_openff.models.OpenMMFF.from_file(amber99sb)
    assert isinstance(ff.data, ff.dtype())
