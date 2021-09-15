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
from mmic_openff.mmic_openff import __version__

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
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        raise NotImplementedError

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

        mmol = inputs.schema_object
        ndim = mmol.ndim

        if ndim != 3:
            raise NotImplementedError(
                "{creator} supports only 3D molecules.".format(**provenance_stamp)
            )  # need to double check this

        mol = OffMolecule()
        mol.name = mmol.name
        natoms = len(mmol.symbols)

        if mmol.atom_labels is not None:
            atom_labels = mmol.atom_labels
        else:
            atom_labels = None

        if mmol.atomic_numbers is not None:
            atomic_numbers = mmol.atomic_numbers
        else:
            raise NotImplementedError(
                "{creator} supports only atomic/molecular systems. Molecule.atomic_numbers must be defined.".format(
                    **provenance_stamp
                )  # need to double check this
            )

        if isinstance(mmol.extras, dict):
            extras = mmol.extas
        else:
            extras = {}

        # not yet supported by MMSchema
        fm_charges = extras.get("formal_charges", [0 for _ in range(natoms)])
        aroma = extras.get("is_aromatic", [False for _ in range(natoms)])

        for index, symb in enumerate(mmol.symbols):
            mol.add_atom(
                atomic_number=mmol.atomic_numbers[index],
                name=atom_labels[index] if atom_labels is not None else symb,
                is_aromatic=aroma[index],
                formal_charge=fm_charges[index],
            )

        if mmol.connectivity is not None:
            for (
                i,
                j,
                order,
            ) in mmol.connectivity:

                mol.add_bond(
                    atom1=i,
                    atom2=j,
                    bond_order=order,
                    is_aromatic=False,
                    stereochemistry=None,
                    fractional_bond_order=None,  # Optional argument; Wiberg (or similar) bond order
                )

        if mmol.geometry is not None:
            try:
                geo_units = getattr(openmm_unit, mmol.geometry_units)
            except AttributeError:
                AttributeError(
                    f"Unit {mmol.geometry_units} not supported by {creator}.".format(
                        **provenance_stamp
                    )
                )

            geo = openmm_unit.Quantity(mmol.geometry.reshape(natoms, ndim), geo_units)
            mol._add_conformer(
                geo.in_units_of(getattr(openmm_unit, openff_units["length"]))
            )

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
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        raise NotImplementedError

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
        if off_mol.n_conformers > 1:
            raise NotImplementedError("Multi-conformers not supported.")

        if off_mol.n_conformers == 1:
            conformer = off_mol.conformers[0]
            geo_units = conformer.unit.get_name()
            geo = conformer._value.flatten()
        else:
            geo = None

        atomic_data = [
            (atom.name, atom.atomic_number, atom.mass._value) for atom in off_mol.atoms
        ]
        names, atomic_nums, masses = zip(*atomic_data)

        masses_units = off_mol.atoms[0].mass.unit.get_name()

        connectivity = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order)
            for bond in off_mol.bonds
        ]

        input_dict = {
            "name": off_mol.name,
            "atomic_numbers": atomic_nums,
            "atom_labels": names,
            "geometry": geo,
            "geometry_units": geo_units,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": masses_units,
        }

        return True, TransOutput(
            proc_input=inputs,
            schema_object=Molecule(**input_dict),
            success=True,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )
