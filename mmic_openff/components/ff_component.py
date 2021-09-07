from mmelemental.models import forcefield
from mmic_translator.models import (
    TransInput,
    TransOutput,
)
from mmic_translator.components import TransComponent
from mmic_openff.mmic_openff import units as openmm_units
from typing import List, Tuple, Optional
from collections.abc import Iterable
from mmelemental.util.units import convert
from simtk.openmm.app import forcefield as openmm_ff
import simtk

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

        empty_atom = parmed.topologyobjects.Atom()
        mmff = inputs.schema_object
        pff = parmed.structure.Structure()

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
        mm_units = forcefield.ForceField.get_units()

        atoms = ff._atomTypes
        templates = ff._templates
        templates_ = [
            {resname: [atom.name for atom in residue.atoms]}
            for resname, residue in templates.items()
        ]

        # names-defs = [(atom.name, symbols[atom.type]) for atom in residue.atoms for residue in templates.items()]

        for gen in ff.getGenerators():
            if isinstance(gen, simtk.openmm.app.forcefield.NonbondedGenerator):
                nonbond_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.HarmonicBondGenerator):
                bond_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.HarmonicAngleGenerator):
                angle_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.PeriodicTorsionGenerator):
                dihedral_ff = gen
            else:
                raise NotImplementedError

        # Need to map these to potential types
        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()

        data = [
            (
                atom.atomClass,
                atom.element.symbol,
                atom.element.atomic_number,
                atom.mass,
            )
            for _, atom in ff._atomTypes.items()
        ]
        types, symbols, atomic_numbers, masses = zip(*data)
        masses_units = atoms["0"].element.mass.unit.get_name()
        masses = convert(masses, masses_units, mm_units["masses_units"])

        nonbonded, charges = self._get_nonbonded(nonbond_ff)
        bonds = self._get_bonds(bond_ff)
        angles = self._get_angles(angle_ff)
        dihedrals = self._get_dihedrals_proper(dihedral_ff)

        # charge_groups = None ... should support charge_groups?
        exclusions = None
        inclusions = None

        input_dict = {
            "masses": masses,
            "charges": charges,
            "bonds": bonds,
            "angles": angles,
            "dihedrals": dihedrals,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "defs": types,  # or names?
            "symbols": symbols,
            # "substructs": residues,
            "atomic_numbers": atomic_numbers,
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

    def _get_nonbonded(self, nonbond):
        lj_scale = nonbond.lj14scale
        coul_scale = nonbond.coulomb14scale

        params = {}
        # How to deal with lj/coul 1-4 scale? relevant to sigma/epsilon scale in parmed?

        for paramName in nonbond.params.paramNames:
            params[paramName] = [
                nonbond.params.paramsForType[str(i)][paramName]
                for i in range(len(nonbond.params.paramsForType))
            ]

        lj = forcefield.nonbonded.potentials.LennardJones(
            sigma=params["sigma"], epsilon=params["epsilon"]
        )

        # Need to include sigma_14 and epsilon_14
        nonbonded = forcefield.nonbonded.NonBonded(params=lj, form="LennardJones")

        # How do we access units?
        return nonbonded, params["charge"]

    def _get_bonds(self, bonds):
        ff = bonds.ff
        bonds_units = forcefield.bonded.bonds.potentials.harmonic.Harmonic.get_units()
        bonds_units.update(forcefield.bonded.Bonds.get_units())

        bonds_lengths = bonds.length
        bonds_k = bonds.k
        ntypes = len(bonds_k)

        connectivity = [
            (
                ff._atomTypes[next(iter(bonds.types1[i]))].atomClass,
                ff._atomTypes[next(iter(bonds.types2[i]))].atomClass,
                1,
            )
            for i in range(ntypes)
        ]

        params = forcefield.bonded.bonds.potentials.Harmonic(
            spring=bonds_k,
            spring_units=f"{openmm_units['energy']} / {openmm_units['length']}**2",
        )

        return forcefield.bonded.Bonds(
            params=params,
            lengths=bonds_lengths,
            lengths_units=openmm_units["length"],
            connectivity=connectivity,
            form="Harmonic",
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
