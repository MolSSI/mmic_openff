"""
mmic_off
Tactic MMIC for OpenFF/MMSchema translation.
"""

# Add imports here
from . import components, models

# Handle versioneer
from ._version import get_versions
from .mmic_off import units

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions


# Need to read these from toolkit.utils.ToolkitWrapper.toolkit_file_read_formats, etc.
molread_ext_maps = {".sdf": "sdf", ".pdb": "pdb"}
molwrite_ext_maps = {".sdf": "sdf", ".pdb": "pdb"}

ffread_ext_maps = {".xml": "xml", ".offxml": "offxml"}
ffwrite_ext_maps = {".xml": "xml", ".offxml": "offxml"}

_classes_map = {
    "Molecule": models.OpenFFMol,
    "ForceField": models.OpenFFSmirnoff,
}
