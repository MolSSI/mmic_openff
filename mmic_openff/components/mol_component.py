from mmelemental.models import Molecule
from openff.toolkit.topology.molecule import Molecule as OffMolecule
from mmic_openff.mmic_openff import units as openff_units
from simtk import (
    unit as openmm_unit,
)  # OpenMM should not be a requirement
from typing import List, Tuple, Optional, Set
from cmselemental.util.decorators import classproperty
from mmic_translator import (
    TransInput,
    TransOutput,
)
from mmic.components import TacticComponent
from mmic_openff.mmic_openff import (
    __version__,
    _supported_versions,
    _mmschema_max_version,
)

provenance_stamp = {
    "creator": "mmic_openff",
    "version": __version__,
    "routine": __name__,
}

__all__ = ["MolToOpenFFComponent", "OpenFFToMolComponent"]


class MolToOpenFFComponent(TacticComponent):
    """A component for converting MMSchema to OpenFF Molecule."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classmethod
    def get_version(cls) -> str:
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

        if inputs.schema_name != Molecule.default_schema_name:
            raise NotImplementedError(
                f"Schema {inputs.schema_name} is not yet supported."
            )

        if inputs.schema_version > _mmschema_max_version:
            raise NotImplementedError(
                f"Schema version {inputs.schema_version} not yet supported."
            )

        mm_mol = inputs.schema_object
        ndim = mm_mol.ndim

        if ndim != 3:
            raise NotImplementedError(
                "{creator} supports only 3D molecules.".format(**provenance_stamp)
            )  # need to double check this

        mol = OffMolecule()
        mol.name = mm_mol.name
        natoms = len(mm_mol.symbols)

        if mm_mol.atomic_numbers is None:
            raise NotImplementedError(
                "{creator} supports only atomic/molecular systems. Molecule.atomic_numbers must be defined.".format(
                    **provenance_stamp
                )  # need to double check this
            )

        extras = getattr(mm_mol, "extras", {}) or {}
        # For now, get any field not supported by MMSchema from extras
        off_mol_unsupport = extras.get(provenance_stamp["creator"], {})
        is_aromatic = off_mol_unsupport.get("is_aromatic")
        stereochem = off_mol_unsupport.get("stereochemistry")
        frac_bond_order = off_mol_unsupport.get("fractional_bond_order")

        if mm_mol.formal_charges is not None:
            formal_charges = mm_mol.formal_charges
            try:
                qunit = getattr(openmm_unit, mm_mol.formal_charges_units)
            except AttributeError:
                raise AttributeError(
                    f"OpenMM unit does not support {mm_mol.formal_charges_units}"
                )
            formal_charges = openmm_unit.Quantity(
                value=formal_charges, unit=mm_mol.qunit
            )
        else:
            formal_charges = None

        for index, symb in enumerate(mm_mol.symbols):
            mol.add_atom(
                atomic_number=int(
                    mm_mol.atomic_numbers[index]
                ),  # must convert numpy.int to int
                name="" if mm_mol.atom_labels is None else mm_mol.atom_labels[index],
                is_aromatic=False if is_aromatic is None else is_aromatic[index],
                formal_charge=0
                if mm_mol.formal_charges is None
                else formal_charges[
                    index
                ],  # what should the formal charge be if unset?
            )

        if mm_mol.connectivity is not None:
            for (
                index,
                (i, j, order),
            ) in enumerate(mm_mol.connectivity):

                mol.add_bond(
                    atom1=i,
                    atom2=j,
                    bond_order=order,
                    is_aromatic=False,  # assume always False for now
                    stereochemistry=None if stereochem is None else stereochem[index],
                    fractional_bond_order=None
                    if frac_bond_order is None
                    else frac_bond_order[index],
                )

        if mm_mol.geometry is not None:
            try:
                geo_units = getattr(openmm_unit, mm_mol.geometry_units)
            except AttributeError:
                AttributeError(
                    f"Unit {mm_mol.geometry_units} not supported by {creator}.".format(
                        **provenance_stamp
                    )
                )

            geo = openmm_unit.Quantity(mm_mol.geometry.reshape(natoms, ndim), geo_units)
            mol._add_conformer(
                geo.in_units_of(getattr(openmm_unit, openff_units["length"]))
            )

        if mm_mol.partial_charges is not None:
            partial_charges = mm_mol.partial_charges
            try:
                qunit = getattr(openmm_unit, mm_mol.partial_charges_units)
            except AttributeError:
                raise AttributeError(
                    f"OpenMM unit does not support {mm_mol.partial_charges_units}"
                )
            partial_charges = openmm_unit.Quantity(value=partial_charges, unit=qunit)
            mol.partial_charges = partial_charges

        success = True
        return success, TransOutput(
            proc_input=inputs,
            data_object=mol,
            success=success,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )


class OpenFFToMolComponent(TacticComponent):
    """A component for converting OpenFF to MMSchema Molecule object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classmethod
    def get_version(cls) -> str:
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

        off_mol = inputs.data_object
        single_atom = off_mol.atoms[0]
        input_dict = {}

        if off_mol.n_conformers > 1:
            raise NotImplementedError("Multi-conformers not supported.")

        if off_mol.n_conformers == 1:
            conformer = off_mol.conformers[0]
            input_dict["geometry"] = conformer._value.flatten()
            input_dict["geometry_units"] = conformer.unit.get_name()

        atomic_data = [
            (atom.name, atom.atomic_number, atom.mass._value) for atom in off_mol.atoms
        ]
        names, atomic_nums, masses = zip(*atomic_data)
        masses_units = single_atom.mass.unit.get_name()  # assume unit does not change

        if single_atom.partial_charge:
            input_dict["partial_charges"] = [
                atom.partial_charge._value for atom in off_mol.atoms
            ]
            input_dict[
                "partial_charges_units"
            ] = single_atom.partial_charge.unit.get_symbol()

        if single_atom.formal_charge:
            input_dict["formal_charges"] = [
                atom.formal_charge._value for atom in off_mol.atoms
            ]
            input_dict[
                "formal_charges_units"
            ] = single_atom.formal_charge.unit.get_symbol()

        connectivity = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order)
            for bond in off_mol.bonds
        ]

        # Deal with data not yet suppported by MMSchema
        bond_extra = [
            (bond.stereochemistry, bond.fractional_bond_order) for bond in off_mol.bonds
        ]
        bonds_stereo, bonds_frac_order = zip(*bond_extra)
        is_aromatic = [atom.is_aromatic for atom in off_mol.atoms]

        input_dict.update(
            {
                "atomic_numbers": atomic_nums,
                "atom_labels": names,
                "connectivity": connectivity,
                "masses": masses,
                "masses_units": masses_units,
                "extras": {
                    provenance_stamp["creator"]: {
                        "stereochemistry": bonds_stereo,
                        "fractional_bond_order": bonds_frac_order,
                        "is_aromatic": is_aromatic,
                    }
                },
            }
        )

        return True, TransOutput(
            proc_input=inputs,
            schema_object=Molecule(**input_dict),
            success=True,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )
