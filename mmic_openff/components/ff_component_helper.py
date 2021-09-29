from mmelemental.models import forcefield
from mmelemental.util import units
from openff.toolkit.utils import string_to_quantity, quantity_to_string
from openmm import unit as openmm_unit
from typing import Dict, Any
from .potential_maps import (
    _dihedrals_potentials_map,
    _dihedrals_improper_potentials_map,
)


def pint_to_openmm(quant: "pint.unit.Quantity") -> openmm_unit.Quantity:
    _, phys_units = quant.to_tuple()
    openmm_quant = openmm_unit.Quantity(quant.magnitude)

    try:
        for unit, exp in phys_units:
            openmm_quant *= getattr(openmm_unit, unit) ** exp

    except AttributeError:
        raise AttributeError(f"OpenMM unit does not support {str(quant.u)}.")

    return openmm_quant


class FFToOpenFFComponentHelper:
    @classmethod
    def _get_bonds(cls, bonds: Dict[str, Any]) -> Dict[str, Any]:
        unsupported_bonds = bonds.extras.get(__package__)

        if unsupported_bonds is None:
            raise NotImplementedError(
                f"{cls.__name__} currently does not work without additional OFF data: id, fractional_bondorder_method, ..."
            )

        springs = [
            quantity_to_string(
                pint_to_openmm(units.Quantity(spring, units=bonds.params.spring_units))
            )
            for spring in bonds.params.spring
        ]

        lengths = [
            quantity_to_string(
                pint_to_openmm(units.Quantity(length, units=bonds.lengths_units))
            )
            for length in bonds.lengths
        ]

        return {
            "Bond": [
                {
                    "smirks": bonds.defs[i],
                    "id": unsupported_bonds.get("id")[i],
                    "k": springs[i],
                    "length": lengths[i],
                }
                for i in range(len(bonds.defs))
            ],
            "version": bonds.version,
            "potential": bonds.form.lower(),  # OFF potentials seem to be always in lower-case
            "fractional_bondorder_method": unsupported_bonds.get(
                "fractional_bondorder_method"
            ),
            "fractional_bondorder_interpolation": unsupported_bonds.get(
                "fractional_bondorder_interpolation"
            ),
        }

    @classmethod
    def _get_angles(cls, angles: Dict[str, Any]) -> Dict[str, Any]:
        unsupported_angles = angles.extras.get(__package__)

        if unsupported_angles is None:
            raise NotImplementedError(
                f"{cls.__name__} currently does not work without additional OFF data: id."
            )

        springs = [
            quantity_to_string(
                pint_to_openmm(units.Quantity(spring, units=angles.params.spring_units))
            )
            for spring in angles.params.spring
        ]

        angles_ = [
            quantity_to_string(
                pint_to_openmm(units.Quantity(angle, units=angles.angles_units))
            )
            for angle in angles.angles
        ]

        return {
            "Angle": [
                {
                    "smirks": angles.defs[i],
                    "angle": angles_[i],
                    "k": springs[i],
                    "id": unsupported_angles.get("id")[i],
                }
                for i in range(len(angles.defs))
            ],
            "version": angles.version,
            "potential": angles.form.lower(),  # OFF potentials seem to be always in lower-case
        }


class OpenFFToFFComponentHelper:
    @classmethod
    def _get_nonbonded(cls, vdW: Dict[str, Any]):

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

    @classmethod
    def _get_bonds(cls, bonds: Dict[str, Any]):
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

    @classmethod
    def _get_angles(cls, angles: Dict[str, Any]):
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

    @classmethod
    def _get_dihedrals_proper(cls, dihedrals: Dict[str, Any]):
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

    @classmethod
    def _get_dihedrals_improper(cls, imdihedrals: Dict[str, Any]):
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
