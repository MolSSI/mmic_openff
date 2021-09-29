from mmelemental.models import forcefield
from mmelemental.util import units
from mmelemental.util.units import convert
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

_dihedrals_potentials_map = {
    "k*(1+cos(periodicity*theta-phase))": "CharmmMulti",
    # need to add all supported potentials in OpenFFTk
}

_dihedrals_improper_potentials_map = {
    "k*(1+cos(periodicity*theta-phase))": "CharmmMulti",
    # need to add all supported potentials in OpenFFTk
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
        try:
            charge_units = getattr(openmm_unit, mmff.charges_units)
        except AttributeError:
            raise AttributeError(f"OpenMM unit does not support {mmff.charges_units}.")

        off_charges = [
            quantity_to_string(openmm_unit.Quantity(charge, unit=charge_units))
            for charge in mmff.charges
        ]
        charges = [
            {"smirks": mmff.defs[i], "charge1": off_charges[i], "name": mmff.symbols[i]}
            for i in range(len(off_charges))
        ]  # need to distinguish name from id

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
                # "Bonds": ...,
                # "Angles": ...,
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
        nonbonded = self._get_nonbonded(vdW_data)
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
        bonds = self._get_bonds(bonds_data)
        angles_data = ff_data["Angles"]
        angles = self._get_angles(angles_data)
        dihedrals_data = ff_data["ProperTorsions"]
        dihedrals = self._get_dihedrals_proper(dihedrals_data)
        dihedrals_improper_data = ff_data["ImproperTorsions"]
        dihedrals_improper = self._get_dihedrals_improper(dihedrals_improper_data)

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

    def _get_nonbonded(self, vdW: Dict[str, Any]):

        # Need to read scale12, scale13, ...
        if vdW["potential"] != "Lennard-Jones-12-6":
            raise NotImplementedError(
                "mmic_openff supports only Lennard-Jones-12-6 potential for now."
            )

        atoms = vdW["Atom"]

        scaling_factor = 2.0 / 2.0 ** (1.0 / 6.0)  # rmin_half = 2^(1/6) * sigma / 2

        data = [
            (
                atom["smirks"],
                atom["id"],
                string_to_quantity(atom["epsilon"])._value,
                string_to_quantity(atom["rmin_half"])._value * scaling_factor
                if "rmin_half" in atom
                else string_to_quantity(atom.get("sigma"))._value,
            )
            for atom in atoms
        ]
        defs, ids, epsilon, sigma = zip(*data)

        # Pint raises UndefinedUnitError if this fails
        single_atom = atoms[0]  # does the unit change with every atom?
        epsilon_units = str(units.Quantity(single_atom["epsilon"]).u)
        sigma_units = str(
            units.Quantity(single_atom.get("rmin_half", single_atom.get("sigma"))).u
        )

        lj = forcefield.nonbonded.potentials.LennardJones(
            sigma=sigma,
            sigma_units=sigma_units,
            epsilon=epsilon,
            epsilon_units=epsilon_units,
        )

        nonbonded = forcefield.nonbonded.NonBonded(
            params=lj,
            form="LennardJones",
            defs=defs,
            combination_rule=vdW["combining_rules"],
            version=vdW["version"],
            extras={
                __package__: {
                    "ids": ids,
                    "method": vdW["method"],
                    "cutoff": vdW["cutoff"],
                    "switch_width": vdW["switch_width"],
                },
            },
        )

        return nonbonded

    def _get_bonds(self, bonds: Dict[str, Any]):
        potential = getattr(
            forcefield.bonded.bonds.potentials, bonds["potential"].capitalize(), None
        )  # MMSchema always uses capitalized potential name
        if not potential:
            raise NotImplementedError(
                f"Potential {bonds['potential']} not supported by MMSchema."
            )

        bonds_units = potential.default_units
        bonds_units.update(forcefield.bonded.Bonds.default_units)

        data = [
            (
                bond["smirks"],
                bond["id"],
                string_to_quantity(bond["length"])._value,
                string_to_quantity(bond["k"])._value,
            )
            for bond in bonds["Bond"]
        ]
        defs, ids, lengths, springs = zip(*data)

        # Pint raises UndefinedUnitError if this fails
        single_bond = bonds["Bond"][0]  # assume the unit remains constant
        spring_units = str(units.Quantity(single_bond["k"]).u)
        lengths_units = str(units.Quantity(single_bond["length"]).u)

        params = potential(spring=springs, spring_units=spring_units)

        return forcefield.bonded.Bonds(
            version=bonds["version"],
            params=params,
            defs=defs,
            lengths=lengths,
            lengths_units=lengths_units,
            form=potential.__name__,
            extras={
                __package__: {
                    "fractional_bondorder_method": bonds["fractional_bondorder_method"],
                    "fractional_bondorder_interpolation": bonds[
                        "fractional_bondorder_interpolation"
                    ],
                    "id": ids,
                },
            },
        )

    def _get_angles(self, angles: Dict[str, Any]):
        potential = getattr(
            forcefield.bonded.angles.potentials, angles["potential"].capitalize(), None
        )  # MMSchema always uses capitalized potential names
        if not potential:
            raise NotImplementedError(
                f"Potential {angles['potential']} not supported by MMSchema."
            )

        # Get units for all physical properties
        angles_units = potential.default_units
        angles_units.update(forcefield.bonded.Angles.default_units)

        data = [
            (
                angle["smirks"],
                angle["id"],
                string_to_quantity(angle["angle"])._value,
                string_to_quantity(angle["k"])._value,
            )
            for angle in angles["Angle"]
        ]
        defs, ids, angles_, springs = zip(*data)

        # Pint raises an UndefinedUnitError if this fails
        single_angle = angles["Angle"][0]  # assume the unit remains constant
        spring_units = str(units.Quantity(single_angle["k"]).u)
        angles_units = str(units.Quantity(single_angle["angle"]).u)

        params = potential(spring=springs, spring_units=spring_units)

        return forcefield.bonded.Angles(
            version=angles["version"],
            params=params,
            defs=defs,
            angles=angles_,
            angles_units=angles_units,
            form=potential.__name__,
            extras={
                __package__: {
                    "id": ids,
                }
            },
        )

    def _get_dihedrals_proper(self, dihedrals: Dict[str, Any]):
        potential_name = _dihedrals_potentials_map.get(dihedrals["potential"], "")
        potential = getattr(
            forcefield.bonded.dihedrals.potentials, potential_name, None
        )
        if not potential:
            raise NotImplementedError(
                f"Potential {dihedrals['potential']} not supported by MMSchema."
            )

        # Get units for all physical properties
        dihedrals_units = potential.default_units
        dihedrals_units.update(forcefield.bonded.Dihedrals.default_units)

        # Read specific FF data
        data = [
            (
                dihedral["smirks"],
                dihedral["id"],
            )
            for dihedral in dihedrals["Proper"]
        ]
        defs, ids = zip(*data)

        no_terms = lambda dih, prop: len(
            [... for key in dih.keys() if key.startswith(prop)]
        )  # hackish?
        data = [
            (
                [dihedral[f"idivf{i+1}"] for i in range(no_terms(dihedral, "idivf"))],
                [
                    dihedral[f"periodicity{i+1}"]
                    for i in range(no_terms(dihedral, "periodicity"))
                ],
                [
                    string_to_quantity(dihedral[f"phase{i+1}"])._value
                    for i in range(no_terms(dihedral, "phase"))
                ],
                [
                    string_to_quantity(dihedral[f"k{i+1}"])._value
                    for i in range(no_terms(dihedral, "k"))
                ],
            )
            for dihedral in dihedrals["Proper"]
        ]
        idivf, periodicity, phase, energy = zip(*data)

        # Pint raises UndefinedUnitError if this fails
        single_dihedral = dihedrals["Proper"][0]  # assume the unit remains constant
        energy_units = str(units.Quantity(single_dihedral["k1"]).u)
        phase_units = str(units.Quantity(single_dihedral["phase1"]).u)

        params = potential(
            energy=energy,
            energy_units=energy_units,
            periodicity=periodicity,
            phase=phase,
            phase_units=phase_units,
        )

        return forcefield.bonded.Dihedrals(
            version=dihedrals["version"],
            params=params,
            defs=defs,
            form=potential.__name__,
            extras={
                __package__: {
                    "fractional_bondorder_method": dihedrals[
                        "fractional_bondorder_method"
                    ],
                    "fractional_bondorder_interpolation": dihedrals[
                        "fractional_bondorder_interpolation"
                    ],
                    "default_idivf": dihedrals["default_idivf"],
                    "idivf": idivf,  # what does this represent?
                },
            },
        )

    def _get_dihedrals_improper(self, imdihedrals: Dict[str, Any]):
        potential_name = _dihedrals_improper_potentials_map.get(
            imdihedrals["potential"], ""
        )
        potential = getattr(
            forcefield.bonded.dihedrals_improper.potentials, potential_name, None
        )
        if not potential:
            raise NotImplementedError(
                f"Potential {dihedrals['potential']} not supported by MMSchema."
            )

        # Get units for all physical properties
        imdihedrals_units = potential.default_units
        imdihedrals_units.update(forcefield.bonded.DihedralsImproper.default_units)

        # Read specific FF data
        data = [
            (
                imdihedral["smirks"],
                imdihedral["id"],
            )
            for imdihedral in imdihedrals["Improper"]
        ]
        defs, ids = zip(*data)

        no_terms = lambda imdih, prop: len(
            [... for key in imdih.keys() if key.startswith(prop)]
        )  # hackish?
        data = [
            (
                [
                    imdihedral[f"periodicity{i+1}"]
                    for i in range(no_terms(imdihedral, "periodicity"))
                ],
                [
                    string_to_quantity(imdihedral[f"phase{i+1}"])._value
                    for i in range(no_terms(imdihedral, "phase"))
                ],
                [
                    string_to_quantity(imdihedral[f"k{i+1}"])._value
                    for i in range(no_terms(imdihedral, "k"))
                ],
            )
            for imdihedral in imdihedrals["Improper"]
        ]
        periodicity, phase, energy = zip(*data)

        # Pint raises UndefinedUnitError if this fails
        single_imdihedral = imdihedrals["Improper"][
            0
        ]  # assume the unit remains constant
        energy_units = str(units.Quantity(single_imdihedral["k1"]).u)
        phase_units = str(units.Quantity(single_imdihedral["phase1"]).u)

        params = potential(
            energy=energy,
            energy_units=energy_units,
            periodicity=periodicity,
            phase=phase,
            phase_units=phase_units,
        )

        return forcefield.bonded.Dihedrals(
            version=imdihedrals["version"],
            params=params,
            defs=defs,
            form=potential.__name__,
            extras={
                __package__: {
                    "default_idivf": imdihedrals["default_idivf"],
                },
            },
        )
