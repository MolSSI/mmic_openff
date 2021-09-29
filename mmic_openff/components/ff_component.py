from mmelemental.models import forcefield
from mmelemental.util import units
from mmic_translator import (
    InputTrans,
    OutputTrans,
)
from mmic.components import TacticComponent
from cmselemental.util.decorators import classproperty
from openmm import unit as openmm_unit
from openff.toolkit.typing.engines import smirnoff
from openff.toolkit.utils import string_to_quantity, quantity_to_string
from openff.toolkit.utils.exceptions import SMIRNOFFParseError
from collections.abc import Iterable
from typing import List, Tuple, Optional, Dict, Any, Set
from .ff_component_helper import (
    OpenFFToFFComponentHelper,
    FFToOpenFFComponentHelper,
    pint_to_openmm,
)

from mmic_openff.mmic_openff import (
    __package__,
    __version__,
    _supported_versions,
    _mmschema_max_version,
)

provenance_stamp = {
    "creator": __package__,
    "version": __version__,
    "routine": __name__,
}

__all__ = ["FFToOpenFFComponent", "OpenFFToFFComponent"]


class FFToOpenFFComponent(TacticComponent):
    """A component for converting MMSchema to OpenFF ForceField object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
        """Returns distutils-style versions of openff.toolkit this component supports.

        Examples
        --------
        The string ">1.0, !=1.5.1, <2.0" implies any version after 1.0 and before 2.0
        is compatible, except 1.5.1

        Returns
        -------
        str
            Return a dist-utils valid version string.

        """
        return _supported_versions

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.

        Returns
        -------
        Set[str]

        """
        return {"mmic_translator"}

    def execute(
        self,
        inputs: InputTrans,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:

        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        mmff = inputs.schema_object
        off = smirnoff.ForceField()
        extras = getattr(mmff, "extras", {}) or {}
        unsupported = extras.get(__package__)

        if unsupported is None:
            raise NotImplementedError(
                f"{self.__repr_name__()} currently does not work without additional OFF data: version, constraints, ..."
            )

        # Get library charges
        off_charges = [
            quantity_to_string(
                pint_to_openmm(units.Quantity(charge, units=mmff.charges_units))
            )
            for charge in mmff.charges
        ]
        charges = [
            {"smirks": mmff.defs[i], "charge1": off_charges[i], "name": mmff.symbols[i]}
            for i in range(len(off_charges))
        ]  # need to distinguish name from id

        bonds = FFToOpenFFComponentHelper._get_bonds(mmff.bonds)
        angles = FFToOpenFFComponentHelper._get_angles(mmff.angles)

        smirnoff_data = {
            "SMIRNOFF": {
                # Metadata
                "version": mmff.version,  # what should be the default version be when unspecified?
                "Date": unsupported.get("Date"),
                "Author": mmff.author,
                # Data
                "aromaticity_model": unsupported.get(
                    "aromaticity_model", "OEAroModel_MDL"
                ),  # must generalize this
                "Constraints": unsupported.get("Constraints"),
                "Bonds": bonds,
                "Angles": angles,
                # "ProperTorsions": ...,
                # "ImproperTorsions": ...,
                # "vdW": ...,
                "Electrostatics": unsupported.get("Electrostatics"),
                "LibraryCharges": {
                    "LibraryCharge": charges,
                    "version": mmff.version,
                },  # use mmff.version for libcharges?
                "ToolkitAM1BCC": unsupported.get("ToolkitAM1BCC"),
            },
        }

        try:
            off._load_smirnoff_data(smirnoff_data)
            success = True
        except SMIRNOFFParseError:
            success = False

        return success, OutputTrans(
            proc_input=inputs,
            data_object=off,
            success=success,
            provenance=provenance_stamp,
            schema_name=mmff.schema_name,
            schema_version=mmff.schema_version,
        )


class OpenFFToFFComponent(TacticComponent):
    """A component for converting OpenFF ForceField to MMSchema object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
        """Returns distutils-style versions of openff.toolkit this component supports.

        Examples
        --------
        The string ">1.0, !=1.5.1, <2.0" implies any version after 1.0 and before 2.0
        is compatible, except 1.5.1

        Returns
        -------
        str
            Return a dist-utils valid version string.

        """
        return _supported_versions

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.

        Returns
        -------
        Set[str]
        """
        return {"mmic_translator"}

    def execute(
        self,
        inputs: InputTrans,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:

        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        ff = inputs.data_object
        ff_data = ff._to_smirnoff_data().get("SMIRNOFF")
        mm_units = forcefield.ForceField.default_units

        # VdW/coulomb data
        vdW_data = ff_data["vdW"]
        nonbonded = OpenFFToFFComponentHelper._get_nonbonded(vdW_data)
        libdata = [
            (
                item["smirks"],
                string_to_quantity(item["charge1"])._value,
                item.get("name", item.get("id")),
            )  # how to distinguish between id and name?
            for item in ff_data["LibraryCharges"]["LibraryCharge"]
        ]
        defs, charges, names = zip(*libdata)
        single_charge = ff_data["LibraryCharges"]["LibraryCharge"][0]
        charges_units = str(units.Quantity(single_charge["charge1"]).u)

        # Bonded data
        bonds_data = ff_data["Bonds"]
        bonds = OpenFFToFFComponentHelper._get_bonds(bonds_data)
        angles_data = ff_data["Angles"]
        angles = OpenFFToFFComponentHelper._get_angles(angles_data)
        dihedrals_data = ff_data["ProperTorsions"]
        dihedrals = OpenFFToFFComponentHelper._get_dihedrals_proper(dihedrals_data)
        dihedrals_improper_data = ff_data["ImproperTorsions"]
        dihedrals_improper = OpenFFToFFComponentHelper._get_dihedrals_improper(
            dihedrals_improper_data
        )

        exclusions = None  # need to read these
        inclusions = None  # not available in OpenFF?

        input_dict = {
            "name": getattr(ff, "name", None),
            "version": ff_data.get(
                "version"
            ),  # not the right version, should be read from the offxml file?
            "author": ff.author,
            "charges": charges,  # librarycharges
            "charges_units": charges_units,
            "defs": defs,
            "symbols": names,
            "bonds": bonds,
            "angles": angles,
            "dihedrals": dihedrals,
            "dihedrals_improper": dihedrals_improper,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "extras": {
                __package__: {
                    "date": ff.date,
                    "aromaticity_model": ff.aromaticity_model,
                    "ToolkitAM1BCC": ff_data[
                        "ToolkitAM1BCC"
                    ],  # for ambertools -> MMSchema agnostic
                    "Electrostatics": ff_data[
                        "Electrostatics"
                    ],  # sim input parameters -> MMSchema agnostic
                    "Constraints": ff_data[
                        "Constraints"
                    ],  # sim input params -> MMSchema agnostic
                },
            },
        }

        ff = forcefield.ForceField(**input_dict)
        success = True
        return success, OutputTrans(
            proc_input=inputs,
            schema_object=ff,
            success=success,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )
