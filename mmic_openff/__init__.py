"""
mmic_openff
Tactic MMIC for OpenFF/MMSchema translation.
"""

# Add imports here
from . import components, models

# Handle versioneer
from ._version import get_versions
from .mmic_openff import units
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions


def _gen_ext_maps(obj, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
    toolkit = None
    ext_map = {}

    for query_toolkit in toolkit_registry.registered_toolkits:
        ext_map.update(
            {"." + val.lower(): val.lower() for val in getattr(query_toolkit, obj)}
        )

    return ext_map


# Need to read these from toolkit.utils.ToolkitWrapper.toolkit_file_read_formats, etc.
molread_ext_maps = _gen_ext_maps("toolkit_file_read_formats")
molwrite_ext_maps = _gen_ext_maps("toolkit_file_write_formats")

ffread_ext_maps = {".xml": "xml", ".offxml": "offxml"}
ffwrite_ext_maps = {".xml": "xml", ".offxml": "offxml"}


_classes_map = {
    "Molecule": models.OpenFFMol,
    "ForceField": models.OpenFFSmirnoff,
}
