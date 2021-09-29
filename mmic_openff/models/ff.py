from typing import Dict, Any, Optional, Union, List
from mmic_translator.models import ToolkitModel
from mmelemental.models import ForceField
from cmselemental.util.decorators import classproperty

# Import Components
from mmic_openff.components.ff_component import FFToOpenFFComponent
from mmic_openff.components.ff_component import OpenFFToFFComponent

# Import OpenFF stuff
from openff.toolkit import __version__ as off_version
from openff.toolkit.typing.engines.smirnoff import ForceField as ForceFieldSmirnoff

__all__ = ["OpenFFSmirnoff"]


class OpenFFSmirnoff(ToolkitModel):
    """A model for storing a Smirnoff ForceField object in the OpenFF toolkit."""

    @classproperty
    def engine(cls):
        return "openff", off_version

    @classproperty
    def dtype(cls):
        """Returns the fundamental force field object type."""
        return ForceFieldSmirnoff

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
        if hasattr(data, "author") and hasattr(data, "date"):
            if callable(data._to_smirnoff_data):  # choose a method specific to SMIRNOFF
                assert data._to_smirnoff_data().get(
                    "SMIRNOFF"
                ), "Could not extract SMIRNOFF data from force field object."
                return data
        raise ValueError("OpenFF Forcefield object is invalid!")

    @classmethod
    def from_file(cls, filename: str, **kwargs) -> "OpenFFSmirnoff":
        """
        Constructs an OpenFF object from file(s).

        Parameters
        ----------
        filename : str
            The forcefield filename to read
        **kwargs
            Any additional keywords to pass to the constructor

        Returns
        -------
        OpenFFSmirnoff
            A constructed OpenFFSmirnoff object.

        """
        ff = cls.dtype(filename, **kwargs)

        return cls(data=ff)

    @classmethod
    def from_schema(
        cls, data: ForceField, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "OpenFFSmirnoff":
        """
        Constructs an OpenFFSmirnoff object from an MMSchema ForceField object.
        Parameters
        ----------
        data: ForceField
            Data to construct the forcefield object from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        OpenMMFF
            A constructed OpenMMFF object.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
            "schema_name": data.schema_name,
            "keywords": kwargs,
        }
        out = FFToOpenFFComponent.compute(inputs)
        return cls(data=out.data_object, data_units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the forcefield to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors.
        """
        if dtype:
            kwargs["format"] = dtype
        self.data.save(filename, **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> ForceField:
        """Converts the forcefield to MMSchema ForceField object.
        Parameters
        ----------
        version: int, optional
            Schema specification version to comply with e.g. 1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {
            "data_object": self.data,
            "schema_version": version,
            "schema_name": kwargs.pop("schema_name", ForceField.default_schema_name),
            "keywords": kwargs,
        }
        out = OpenMMToFFComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
