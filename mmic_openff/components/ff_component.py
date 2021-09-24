from mmelemental.models import forcefield
from mmelemental.util import units
from mmelemental.util.units import convert
from mmic_translator import (
    TransInput,
    TransOutput,
)
from mmic.components import TacticComponent
from cmselemental.util.decorators import classproperty
from openff.toolkit.typing.engines import smirnoff
from openff.toolkit.utils import string_to_quantity
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

__all__ = ["FFToOpenFFComponent", "OpenFFToFFComponent"]


class FFToOpenFFComponent(TacticComponent):
    """A component for converting MMSchema to OpenFF ForceField object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classproperty
    def version(cls) -> str:
        """Returns distutils-style version string.

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
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        mmff = inputs.schema_object
        off = openff.toolkit.typing.engines.smirnoff.ForceField()
        extras = getattr(mmff, "extras", {}) or {}
        unsupported = extras.get(__package__)

        # Get library charges
        libcharges = mmff.charges

        smirnoff_data = {
            "SMIRNOFF": {
                # Metadata
                "version": mmff.version,
                "Author": mmff.author,
                "Date": unsupported.get("date"),
                "aromaticity_model": unsupported.get("aromaticity_model"),
                # Data
                "Constraints": unsupported.get("Constraints"),
                "Bonds": ...,
                "Angles": ...,
                "ProperTorsions": ...,
                "ImproperTorsions": ...,
                "vdW": ...,
                "Electrostatics": ...,
                "LibraryCharges": {
                    "LibraryCharge": charges,
                    "version": mmff.version,
                },  # use ff.version for libcharges?
                "ToolkitAM1BCC": unsupported.get(""),
            },
        }

        try:
            off._load_smirnoff_data(smirnoff_data)
        except SMIRNOFFParseError:
            raise SMIRNOFFParseError

        success = True
        return success, TransOutput(
            proc_input=inputs,
            data_object=pff,
            success=success,
            provenance=provenance_stamp,
            schema_name=mmff.schema_name,
            schema_version=mmff.schema_version,
        )


class OpenFFToFFComponent(TacticComponent):
    """A component for converting OpenFF ForceField to MMSchema object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classproperty
    def version(cls) -> str:
        """Returns distutils-style version string.

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
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        ff = inputs.data_object
        ff_data = ff._to_smirnoff_data().get("SMIRNOFF")
        mm_units = forcefield.ForceField.default_units

        vdW_data = ff_data["vdW"]
        bonds_data = ff_data["Bonds"]
        angles_data = ff_data["Angles"]

        nonbonded = self._get_nonbonded(vdW_data)
        bonds = self._get_bonds(bonds_data)
        angles = self._get_angles(angles_data)
        # dihedrals = self._get_dihedrals_proper(dihedral_ff)

        # charge_groups = None ... should support charge_groups?
        exclusions = None
        inclusions = None

        input_dict = {
            "name": getattr(ff, "name", None),
            "version": ff_data["version"],  # not sure abt this
            "author": ff.author,
            # "masses": masses,
            # "charges": charges,
            # "bonds": bonds,
            # "angles": angles,
            # "dihedrals": dihedrals,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            # "defs": types,  # or names?
            # "symbols": symbols,
            # "atomic_numbers": atomic_numbers,
            "extras": {
                __package__: {
                    "date": ff.date,
                    "aromaticity_model": ff.aromaticity_model,
                    "ToolkitAM1BCC": ff_data["ToolkitAM1BCC"],
                },
            },
        }

        ff = forcefield.ForceField(**input_dict)
        success = True
        return success, TransOutput(
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

        # Pint raises an UndefinedUnitError if this fails
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
        try:
            potential = getattr(
                forcefield.bonded.bonds.potentials, bonds["potential"].capitalize()
            )
        except AttributeError:
            raise AttributeError(
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

        # Pint raises an UndefinedUnitError if this fails
        single_bond = bonds["Bond"][0]  # assume the unit remains constant
        spring_units = str(units.Quantity(single_bond["k"]).u)
        lengths_units = str(units.Quantity(single_bond["length"]).u)

        params = potential(spring=springs, spring_units=spring_units)

        return forcefield.bonded.Bonds(
            version=bonds["version"],
            params=params,
            lengths=lengths,
            lengths_units=lengths_units,
            form=potential.__name__,
            extras={
                __package__: {
                    "fractional_bondorder_method": bonds["fractional_bondorder_method"],
                    "fractional_bondorder_interpolation": bonds[
                        "fractional_bondorder_interpolation"
                    ],
                },
            },
        )

    def _get_angles(self, angles: Dict[str, Any]):
        try:
            potential = getattr(
                forcefield.bonded.angles.potentials, angles["potential"].capitalize()
            )
        except AttributeError:
            raise AttributeError(
                f"Potential {angles['potential']} not supported by MMSchema."
            )

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
            angles=angles_,
            angles_units=angles_units,
            form=potential.__name__,
        )

    def _get_dihedrals_proper(self, dihedrals):
        ff = dihedrals.ff
        proper = dihedrals.proper
        dihedrals_units = (
            forcefield.bonded.dihedrals.potentials.harmonic.Harmonic.default_units
        )
        dihedrals_units.update(forcefield.bonded.Dihedrals.default_units)

        connectivity = [
            (
                ff._atomTypes[next(iter(dihedral.types1))].atomClass,
                ff._atomTypes[next(iter(dihedral.types2))].atomClass,
                ff._atomTypes[next(iter(dihedral.types3))].atomClass,
                ff._atomTypes[next(iter(dihedral.types4))].atomClass,
            )
            for dihedral in proper
        ]

        fields = [
            (dihedral.k, dihedral.periodicity, dihedral.phase) for dihedral in proper
        ]
        energy, periodicity, phase = zip(*fields)

        params = forcefield.bonded.dihedrals.potentials.CharmmMulti(
            energy=energy,
            energy_units=openmm_units["energy"],
            periodicity=periodicity,
            phase=phase,
            phase_units=openmm_units["angle"],
        )

        return forcefield.bonded.Dihedrals(
            params=params,
            connectivity=connectivity,
            form="CharmmMulti",
        )
