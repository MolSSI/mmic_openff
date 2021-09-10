from typing import Dict, Any, Optional
from mmic_translator.models.base import ToolkitModel
from mmelemental.models import Molecule
from cmselemental.types import Array
import numpy
from pydantic import Field
from pathlib import Path


# Import OpenFF stuff
from openff.toolkit import __version__ as off_version
from openff.toolkit import topology as off_top

__all__ = ["OpenFFMol"]


class OpenFFMol(ToolkitModel):
    """A model for OpenFF molecule class storing a molecule object."""

    @classmethod
    def engine(cls):
        return "openff", off_version

    @classmethod
    def dtype(cls):
        """Returns the fundamental molecule object type."""
        return off_top.Molecule

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
        if len(data.atoms) > 0:
            return data
        raise ValueError("OpenFF Molecule object does not contain any atoms!")

    @classmethod
    def from_file(
        cls, filename: str, top_filename: str = None, **kwargs
    ) -> "OpenFFMol":
        """
        Constructs an OpenFFMol object from file(s).

        Parameters
        ----------
        filename : str
            The atomic positions filename to read
        top_filename: str, optional
            The topology filename to read
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        OpenFFMol
            A constructed OpenFFMol class.
        """

        if top_filename:
            raise NotImplementedError(
                f"{cls} does not support passing topology filenames"
            )

        mol = cls.dtype().from_file(filename, **kwargs)

        return cls(
            data=mol,
        )

    @classmethod
    def from_schema(
        cls, data: Molecule, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "OpenFFMol":
        """
        Constructs an OpenFFMol object from an MMSchema Molecule object.
        Parameters
        ----------
        data: Molecule
            Data to construct Molecule from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        OpenFFMol
            A constructed OpenFFMol class.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
        }
        out = MolToOpenFFComponent.compute(inputs)
        return cls(data=out.data_object, units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, mode: str = None, **kwargs):
        """Writes the molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors. kwargs takes precedence over  data.
        """
        if mode:
            raise NotImplementedError(f"{cls} does not support specifying write mode.")

        if not dtype:
            dtype = Path(filename).suffix.removeprefix(".")

        self.data.to_file(filename, file_format=dtype, **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> Molecule:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        kwargs["positions"] = kwargs.get("positions", self.positions)
        kwargs["positions_units"] = kwargs.get("positions_units", self.positions_units)

        inputs = {"data_object": self.data, "schema_version": version, "kwargs": kwargs}
        out = OpenFFToMolComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
