"""
mmic_openff
Tactic MMIC for OpenFF/MMSchema translation.
"""

# Add imports here
from . import components, models

# Handle versioneer
from .mmic_openff import __version__, units, _gen_ext_maps

# Need to read these from toolkit.utils.ToolkitWrapper.toolkit_file_read_formats, etc.
molread_ext_maps = _gen_ext_maps("toolkit_file_read_formats")
molwrite_ext_maps = _gen_ext_maps("toolkit_file_write_formats")

ffread_ext_maps = {".xml": "xml", ".offxml": "offxml"}
ffwrite_ext_maps = {".xml": "xml", ".offxml": "offxml"}


_classes_map = {
    "Molecule": models.OpenFFMol,
    "ForceField": models.OpenFFSmirnoff,
}
