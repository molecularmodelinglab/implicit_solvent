import sys

from tqdm import trange
from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
import openmm as mm
import numpy as np
from openmm import app
from openmm import unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule
from openff.units.openmm import to_openmm

# FYI: to fix the PDB file, I used pdbfixer:
# pdbfixer data/5vsc_A_rec.pdb --output data/5vsc_A_rec_fixed.pdb --add-residues --add-atoms=all --keep-heterogens=none


def get_minimized_energy(pdb, lig):

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.add(lig.to_topology().to_openmm(), to_openmm(lig.conformers[0]))
    mergedTopology = modeller.topology
    mergedPositions = modeller.positions

    # load the AMBER force field with implicit solvent
    forcefield_kwargs = {
        "nonbondedCutoff": 2*unit.nanometer,
        "constraints": "HBonds"
    }
    ffs = ['amber/ff14SB.xml', 'implicit/gbn2.xml']
    system_generator = SystemGenerator(forcefields=ffs, small_molecule_forcefield='gaff-2.11', molecules=[lig], forcefield_kwargs=forcefield_kwargs, cache='data/db.json')
    system = system_generator.create_system(mergedTopology)

    residues = list(modeller.topology.residues())
    # freeze alpha carbons (takes a while if we don't and I don't think it matters)
    for r in residues:
        for atom in r.atoms():
            if atom.name == "CA":
                system.setParticleMass(atom.index, 0.0)

    lig_indexes = []
    for atom in residues[-1].atoms():
        lig_indexes.append(atom.index)
    lig_indexes = np.array(lig_indexes)

    # Set up the integrator for energy minimization
    integrator = mm.VerletIntegrator(0.001*unit.picoseconds)

    # Set up the simulation context with the system, integrator, and positions from the PDB file
    context = mm.Context(system, integrator)
    context.setPositions(mergedPositions)

    # Minimize the energy of the system

    print(f'Minimizing energy...')
    mm.LocalEnergyMinimizer.minimize(context)

    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()

    return energy

if __name__ == "__main__":

    crystal_lig_file = "data/4nvq_2od_lig.sdf"
    docked_lig_file = "data/docked.sdf"
    rec_file = "data/5vsc_A_rec_fixed.pdb"

    pdb_openmm = app.PDBFile(rec_file)
    docked_ligs_openmm = Molecule.from_file(docked_lig_file)
    
    # you can get the minimized energy for each docked ligand
    print("Energy of pose 0", get_minimized_energy(pdb_openmm, docked_ligs_openmm[0]))

    docked_ligs = [ mol for mol in Chem.SDMolSupplier(docked_lig_file) ]
    crystal_lig = Chem.SDMolSupplier(crystal_lig_file)[0]

    # you can use the CalcRMS function from rdkit to calculate the RMSD between two molecules
    print("RMSD of pose 0", CalcRMS(crystal_lig, docked_ligs[0]))