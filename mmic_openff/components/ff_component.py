from mmelemental.models import forcefield
from mmelemental.util import units
from mmic_translator import (
    TransComponent,
    TransInput,
    TransOutput,
)
from mmic_openff.mmic_openff import __version__
from typing import List, Tuple, Optional, Dict, Any
from collections.abc import Iterable
from mmelemental.util.units import convert
from openff.toolkit.typing.engines import smirnoff

provenance_stamp = {
    "creator": "mmic_openff",
    "version": __version__,
    "routine": __name__,
}

__all__ = ["FFToOpenFFComponent", "OpenFFToFFComponent"]


class FFToOpenFFComponent(TransComponent):
    """A component for converting MMSchema to OpenFF ForceField object."""

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

        masses = convert(
            mmff.masses, mmff.masses_units, empty_atom.umass.unit.get_symbol()
        )

        charges = TransComponent.get(mmff, "charges")
        charges = convert(
            charges, mmff.charges_units, empty_atom.ucharge.unit.get_symbol()
        )

        atomic_numbers = TransComponent.get(mmff, "atomic_numbers")
        atom_types = TransComponent.get(mmff, "defs")

        rmin, epsilon = self._get_nonbonded(mmff, empty_atom)

        for index, symb in enumerate(mmff.symbols):

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            if atomic_numbers is not None:
                atomic_number = atomic_numbers[index]
            else:
                atomic_number = None

            if atom_types is not None:
                atom_type = atom_types[index]
            else:
                atom_type = None

            if masses is not None:
                mass = masses[index]
            else:
                mass = None

            if charges is not None:
                charge = charges[index]
            else:
                charge = None

            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=symb,
                type=atom_type,
                mass=mass,
                charge=charge,
                nb_idx=0,
                solvent_radius=0.0,
                screen=0.0,
                tree="BLA",
                join=0.0,
                irotat=0.0,
                occupancy=1.0,
                bfactor=0.0,
                altloc="",
                number=-1,
                rmin=rmin[index],
                epsilon=epsilon[index],
                rmin14=None,
                epsilon14=None,
                # bonds=..., faster than connecting atoms one by one as done below?
                # angles=...,
                # dihedrals=...,
                # impropers=...,
                # polarizable=...,
            )

            residues = TransComponent.get(mmff, "substructs")

            if residues:
                resname, resnum = residues[index]
            else:
                raise NotImplementedError(
                    "Residues must be supplied for forcefields based on atom typing."
                )

            pff.add_atom(atom, resname, resnum, chain="", inscode="", segid="")

        # Bonds
        bonds = TransComponent.get(mmff, "bonds")
        if bonds is not None:
            assert (
                mmff.bonds.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                bonds.params.spring, bonds.params.spring_units, "kcal/mol/angstroms**2"
            )
            req = convert(bonds.lengths, bonds.lengths_units, "angstroms")

            for (
                bi,
                (
                    i,
                    j,
                    order,
                ),
            ) in enumerate(mmff.bonds.indices):
                btype = parmed.topologyobjects.BondType(
                    k=spring[bi], req=req[bi], list=pff.bond_types
                )
                pff.bonds.append(
                    parmed.topologyobjects.Bond(
                        pff.atoms[i], pff.atoms[j], order=order, type=btype
                    )
                )
                pff.bond_types.append(btype)
                # both implementations seem to perform almost the same:
                # pff.atoms[i].bond_to(pff.atoms[j])

        # Angles
        angles = TransComponent.get(mmff, "angles")
        if angles is not None:
            assert (
                mmff.angles.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                angles.params.spring, angles.params.spring_units, "kcal/mol/radians^2"
            )
            angles_eq = convert(angles.angles, angles.angles_units, "degrees")

            for ai, (i, j, k) in enumerate(mmff.angles.indices):
                atype = parmed.topologyobjects.AngleType(
                    k=spring[ai], theteq=angles_eq[ai], list=pff.angle_types
                )
                pff.angles.append(
                    parmed.topologyobjects.Angle(
                        pff.atoms[i], pff.atoms[j], pff.atoms[k], type=atype
                    )
                )
                pff.angle_types.append(atype)

        # Dihedrals
        dihedrals = TransComponent.get(mmff, "dihedrals")
        if dihedrals is not None:
            dihedrals = (
                dihedrals.pop() if isinstance(dihedrals, list) else dihedrals
            )  # For now, keep track of only a single type
            # Need to change this ASAP! Must take multiple types into account!
            assert (
                dihedrals.form == "Charmm" or dihedrals.form == "CharmmMulti"
            ), "Only Charmm-style potentials supported for now"

            energy = convert(
                dihedrals.params.energy, dihedrals.params.energy_units, "kcal/mol"
            )
            phase = convert(
                dihedrals.params.phase, dihedrals.params.phase_units, "degrees"
            )
            periodicity = dihedrals.params.periodicity

            for di, (i, j, k, l) in enumerate(dihedrals.indices):
                if isinstance(energy[di], Iterable):
                    dtype = [
                        parmed.topologyobjects.DihedralType(
                            phi_k=energy[di][dj],
                            per=periodicity[di][dj],
                            phase=phase[di][dj],
                            # scee,
                            # scnb,
                            list=pff.dihedral_types,
                        )
                        for dj in range(len(energy[di]))
                    ]
                else:
                    dtype = parmed.topologyobjects.DihedralType(
                        phi_k=energy[di],
                        per=periodicity[di],
                        phase=phase[di],
                        # scee
                        # scnb
                        list=pff.dihedral_types,
                    )
                # assert:
                # dtype.funct = (
                #    9  # hackish: assume all dihedrals are proper and charmm-style
                # )
                pff.dihedrals.append(
                    parmed.topologyobjects.Dihedral(
                        pff.atoms[i],
                        pff.atoms[j],
                        pff.atoms[k],
                        pff.atoms[l],
                        improper=False,
                        type=dtype,
                    )
                )
                pff.dihedral_types.append(dtype)

        smirnoff_data = {
            "SMIRNOFF": {
                # Metadata
                "version": ...,
                "Author": mmff.author,
                "Date": mmff.extras.get("date")
                if isinstance(mmff.extras, dict)
                else None,
                "aromaticity_model": ...,
                # Data fields
                "Constraints": ...,
                "Bonds": ...,
                "Angles": ...,
                "ProperTorsions": ...,
                "ImproperTorsions": ...,
                "vdW": ...,
                "Electrostatics": ...,
                "LibraryCharges": ...,
                "ToolkitAM1BCC": ...,
            },
        }

        return True, TransOutput(proc_input=inputs, data_object=pff)

    def _get_nonbonded(
        self,
        mmff: forcefield.ForceField,
        empty_atom: "parmed.topologyobjects.Atom",
    ) -> Tuple["numpy.ndarray", "numpy.ndarray"]:

        assert (
            mmff.nonbonded.form == "LennardJones"
        ), "Only LJ potential supported for now"

        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()
        scaling_factor = 2 ** (1.0 / 6.0)  # rmin = 2^(1/6) sigma

        rmin = mmff.nonbonded.params.sigma * scaling_factor
        rmin = convert(rmin, lj_units["sigma_units"], empty_atom.urmin.unit.get_name())
        # atom.rmin_14 * rmin_14_factor * scaling_factor,
        epsilon = convert(
            mmff.nonbonded.params.epsilon,
            lj_units["epsilon_units"],
            empty_atom.uepsilon.unit.get_name(),
        )
        # atom.epsilon_14 * epsilon_14_factor,
        return rmin, epsilon


class OpenFFToFFComponent(TransComponent):
    """A component for converting OpenFF ForceField to MMSchema object."""

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
        mm_units = forcefield.ForceField.get_units()

        vdW_data = ff_data["vdW"]
        bonds_data = ff_data["Bonds"]

        nonbonded = self._get_nonbonded(vdW_data)
        # bonds = self._get_bonds(bond_ff)
        # angles = self._get_angles(angle_ff)
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
                provenance_stamp["creator"]: {
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

        scaling_factor = 2.0 / 2.0 ** (1.0 / 6.0)  # rmin_half = 2^(1/6) sigma / 2

        data = [
            (
                atom["smirks"],
                atom["id"],
                float(
                    atom["epsilon"].split("*")[0]
                ),  # hackish but faster than using pint.Quantity
                float(atom.get("rmin_half", atom.get("sigma")).split("*")[0])
                * scaling_factor,  # hackish but faster than using pint.Quantity
            )
            for atom in atoms
        ]
        defs, ids, epsilon, sigma = zip(*data)

        # Pint raises an UndefinedUnitError if this fails
        atom = atoms[0]  # does the unit change with every atom?
        epsilon_units = str(units.Quantity(atom["epsilon"]).u)
        sigma_units = str(units.Quantity(atom.get("rmin_half", atom.get("sigma"))).u)

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
                provenance_stamp["creator"]: {
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

        bonds_units = potential.get_units()
        bonds_units.update(forcefield.bonded.Bonds.get_units())

        data = [
            (
                bond["smirks"],
                bond["id"],
                float(bond["length"].split("*")[0]),
                float(bond["k"].split("*")[0]),
            )
            for bond in bonds["Bond"]
        ]
        defs, ids, lengths, springs = zip(*data)

        # Pint raises an UndefinedUnitError if this fails
        spring_units = str(units.Quantity(bond["length"]).u)
        lengths_units = str(units.Quantity(bond["k"]).u)

        nbonds = len(lengths)
        params = potential(spring=springs, spring_units=spring_units)

        return forcefield.bonded.Bonds(
            version=bonds["version"],
            params=params,
            lengths=bonds_lengths,
            lengths_units=lengths_units,
            form=potential.__name__,
            extras={
                provenance_stamp["creator"]: {
                    "fractional_bondorder_method": bonds["fractional_bondorder_method"],
                    "fractional_bondorder_interpolation": bonds[
                        "fractional_bondorder_interpolation"
                    ],
                },
            },
        )

    def _get_angles(self, angles):
        ff = angles.ff
        angles_units = forcefield.bonded.angles.potentials.harmonic.Harmonic.get_units()
        angles_units.update(forcefield.bonded.Angles.get_units())

        angles_lengths = angles.angle
        angles_k = angles.k
        ntypes = len(angles_k)

        connectivity = [
            (
                ff._atomTypes[next(iter(angles.types1[i]))].atomClass,
                ff._atomTypes[next(iter(angles.types2[i]))].atomClass,
                ff._atomTypes[next(iter(angles.types3[i]))].atomClass,
            )
            for i in range(ntypes)
        ]

        params = forcefield.bonded.angles.potentials.Harmonic(
            spring=angles_k,
            spring_units=f"{openmm_units['energy']} / {openmm_units['angle']}**2",
        )

        return forcefield.bonded.Angles(
            params=params,
            angles=angles_lengths,
            angles_units=openmm_units["angle"],
            connectivity=connectivity,
            form="Harmonic",
        )

    def _get_dihedrals_proper(self, dihedrals):
        ff = dihedrals.ff
        proper = dihedrals.proper
        dihedrals_units = (
            forcefield.bonded.dihedrals.potentials.harmonic.Harmonic.get_units()
        )
        dihedrals_units.update(forcefield.bonded.Dihedrals.get_units())

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
