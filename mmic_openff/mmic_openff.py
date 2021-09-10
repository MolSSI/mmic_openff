"""
mmic_openff.py
A short description of the project.

Handles the primary functions
"""

from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# Incorrect ~ must update
units = {
    "length": "angstrom",
    "time": "ps",
    "energy": "kJ/mol",
    "charge": "e",
    "mass": "amu",
    "angle": "radians",
    "temperature": "kelvin",
}


def _gen_ext_maps(obj, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
    toolkit = None
    ext_map = {}

    for query_toolkit in toolkit_registry.registered_toolkits:
        ext_map.update(
            {"." + val.lower(): val.lower() for val in getattr(query_toolkit, obj)}
        )

    return ext_map
