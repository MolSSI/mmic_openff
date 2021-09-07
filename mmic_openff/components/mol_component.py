from mmelemental.models import Molecule
from openff.toolkit.topology.molecule import Molecule as OffMolecule
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
from mmic_translator import (
    TransComponent,
    TransInput,
    TransOutput,
    __version__,
)
import numpy
from simtk import unit  # Importing OpenMM just to pass geometry units ... how stupid

provenance_stamp = {
    "creator": "mmic_openff",
    "version": __version__,
    "routine": __name__,
}

__all__ = ["MolToOpenFFComponent", "OpenFFToMolComponent"]


class MolToOpenFFComponent(TransComponent):
    """A component for converting MMSchema to OpenFF Molecule."""

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
            raise NotImplementedError("mmic_openff supports only 3D molecules.")

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
                "mmic_openff is supported only for atomic/molecular systems. Molecule.atomic_numbers must be defined."
            )

        if isinstance(mmol.extras, dict):
            extras = mmol.extas
        else:
            extras = {}

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
            geo = unit.Quantity(
                mmol.geometry.reshape(natoms, ndim), getattr(unit, mmol.geometry_units)
            )
            mol._add_conformer(geo.in_units_of(unit.angstrom))

        success = True
        return success, TransOutput(
            proc_input=inputs,
            data_object=mol,
            success=success,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )


class OpenFFToMolComponent(TransComponent):
    """A component for converting OpenFF to MMSchema Molecule object."""

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

        top = inputs.data_object
        geo = inputs.keywords.get("positions", None)
        if geo is not None:
            geo_units = inputs.keywords.get("positions_units", geo.unit.get_name())

            geo = numpy.array([(pos.x, pos.y, pos.z) for pos in geo]).T.flatten()

        atomic_data = [
            (atom.name, atom.element.atomic_number, atom.element.symbol)
            for atom in top.atoms()
        ]
        names, atomic_nums, element_names = zip(*atomic_data)

        try:
            masses = [atom.element.mass._value for atom in top.atoms()]
            masses_units = next(top.atoms()).element.mass.unit.get_name()
        except Exception:
            masses = None
            masses_units = "dalton"

        # If bond order is none, set it to 1.
        connectivity = [
            (bond.atom1.index, bond.atom2.index, bond.order or 1)
            for bond in top.bonds()
        ]

        residues = [(atom.residue.name, atom.residue.index) for atom in top.atoms()]

        input_dict = {
            "atomic_numbers": atomic_nums,
            "symbols": element_names,
            "atom_labels": names,
            "geometry": geo,
            "geometry_units": geo_units,
            "substructs": residues,
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
