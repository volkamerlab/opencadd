"""

This script runs a molecular dynamics simulation for a protein.
In total there are three simulations, starting with an NVT for 10 ns,
then an NPT fpr 10 ns and finally another NVT for 100 ns.
The first NVT is for equilibrating, the NPT for adjusting the pressure.
With the second NVT the final trajectory is built.
Depending on whether the script is run on a cluster or local, the productionSteps must be adapted,
because of high processing demand.
(productionSteps = 10 000 000 fs ('10ns') / 2 fs/step = 5 000 000 step)
The result files are written into the same directory as the input PDB file.

Usage:
python run_md_simulation.py <path to fixed input PDB>

Example:
python run_md_simulation.py ../5IF1/5if1_fixed.pdb

Reference:
Script is adapted from:
http://openmm.org/tutorials/hkmt_tip4pew/

"""
from __future__ import print_function

import os
import sys

import openmm as mm
from openmm import app
from simtk import unit

from sys import stdout

def simulate(
    input_pdb_path,
    output_pdb_path,
    output_trajectory_path,
    production_steps,
    initialize=False,
    checkpoint_in=None,
    checkpoint_out=None,
    extra_force=None,
    use_cluster=True,
):
    """
    Runs a molecular dynamics simulation for a protein.

    Parameters
    ----------
    input_pdb_path : str
        Path to input fixed PDB.
    output_pdb_path : str
        Path to output for PDB files.
    output_trajectory_path : str
        Path to output for trajectories.
    production_steps : int
        Number of simulation steps.
    initialize : bool
        Whether to minimize energies and equilibrate.
    checkpoint_in : str
        Optional starting point for simulation.
    checkpoint_out : str
        Optionally write end state of simulation as checkpoint.
    extra_force : MonteCarloBarostat
        Necessary for NPT.
    use_cluster : bool
        Whether to assume cluster environment or local machine.
        GPU is used on cluster, CPU otherwise.

    """

    # load in input PDB file and force field XML files
    pdb = app.PDBFile(input_pdb_path)
    force_field = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # use app.Modeller to add hydrogen and solvent
    modeller = app.Modeller(pdb.topology, pdb.positions)

    if initialize:
        modeller.addHydrogens(force_field)
        modeller.addSolvent(
            force_field,
            model="tip3p",
            padding=1.0 * unit.nanometers,
            positiveIon="Na+",
            negativeIon="Cl-",
            neutralize=True,
        )

    # prepare system (values from API reference)
    system = force_field.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    # set integrator which is used to simulate the movement of molecules in the system
    integrator = mm.LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
    )
    integrator.setConstraintTolerance(0.00001)

    # add extra force for the MonteCarloBarostat (2. simulation NPT)
    if extra_force:
        system.addForce(extra_force)

    if use_cluster:
        # prepare platform for cluster
        platform = mm.Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision": "mixed"}
    else:
        # prepare platform for laptop
        platform = mm.Platform.getPlatformByName("CPU")
        properties = {}

    # prepare simulation
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    if checkpoint_in:
        with open(checkpoint_in, "rb") as f:
            simulation.context.loadCheckpoint(f.read())

    if initialize:
        # minimize
        print("Minimizing...")
        simulation.minimizeEnergy(maxIterations=1000)

        # equilibrate for 100 steps
        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
        equilibrate_step = max(100, int(production_steps / 1000))
        print("Equilibrating...")
        simulation.step(equilibrate_step)

    # append reporters
    # set how many frames will be written to file
    simulation.reporters.append(app.DCDReporter(output_trajectory_path, 5000))

    # every 50 steps the progress is documented and reported
    simulation.reporters.append(
        app.StateDataReporter(
            stdout,
            50,
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=production_steps,
            separator="\t",
        )
    )
    # write output pdb file and update it every 100 000 steps
    simulation.reporters.append(app.PDBReporter(output_pdb_path, 100000))

    # run simulation with pre-defined settings
    print("Running Production...")
    simulation.step(production_steps)

    if checkpoint_out:
        with open(checkpoint_out, "wb") as f:
            f.write(simulation.context.createCheckpoint())


def run(pdb_filepath):
    output_path_prefix = os.path.splitext(pdb_filepath)[0]

    # 1. NVT simulation 10 ns
    output_path_nvt1_pdb = output_path_prefix + "_nvt1.pdb"
    simulate(
        pdb_filepath,
        output_path_nvt1_pdb,
        output_path_prefix + "_nvt1_trajectory.dcd",
        5000000,
        initialize=True,
        checkpoint_out=output_path_prefix + "_checkpoint_nvt1.chk"
    )
    print("First simulation done")

    # 2. NPT simulation 10 ns
    barostat = mm.MonteCarloBarostat(1 * unit.atmospheres, 300 * unit.kelvin, 50)
    output_path_npt_pdb = output_path_prefix + "_npt.pdb"
    simulate(
        output_path_nvt1_pdb,
        output_path_npt_pdb,
        output_path_prefix + "_npt_trajectory.dcd",
        5000000,
        checkpoint_in=output_path_prefix + "_checkpoint_nvt1.chk",
        checkpoint_out=output_path_prefix + "_checkpoint_npt.chk",
        extra_force=barostat
    )
    print("Second simulation done")

    # 3. NVT simulation 100ns
    simulate(
        output_path_npt_pdb,
        output_path_prefix + "_nvt2.pdb",
        output_path_prefix + "_nvt2_trajectory.dcd",
        50000000,
        checkpoint_in=output_path_prefix + "_checkpoint_npt.chk",
        checkpoint_out=output_path_prefix + "_checkpoint_nvt2.chk"
    )
    print("Third simulation done; Program finished")
    return
