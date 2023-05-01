import openmm as mm
from openmm import app as mmapp
import sys

# Input PDB must contain all the atoms needed by the force field.
pdb = mmapp.PDBFile("/Users/home/Downloads/output.pdb")
# Alternatively, use mmap.PDBxFile

forcefield = mmapp.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=mmapp.PME,
    nonbondedCutoff=1*mm.unit.nanometer,
    constraints=mmapp.HBonds,
)
integrator = mm.LangevinMiddleIntegrator(
    300*mm.unit.kelvin,
    1/mm.unit.picosecond,
    0.004*mm.unit.picoseconds
)

simulation = mmapp.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(mmapp.PDBReporter('output.pdb', 1000))
simulation.reporters.append(
    mmapp.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True)
)
simulation.step(100)
