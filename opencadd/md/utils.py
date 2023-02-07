import os
import parmed

import simtk.openmm as mm
from simtk.openmm import app

from simtk import unit
from pathlib import Path
import pytraj as pt
from pytraj.cluster import kmeans

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

from opencadd._typing import PathLike



def create_topology(filepath_pdf):
    """
    This script creates the topology file needed for clustering with pytraj. It is recommended to use the nvt1 or npt PDB
    file from the previous MD simulation for this step, because these are already reduced in size and this step will take
    less time.

    The topology file defines which atoms are connected to one another through chemical bonds. It  specify
    bonds (2 atoms connected), angles (3 atoms connected) and dihedrals (4 atoms connected linearly).

    Usage:
    python create_topology.py <path to input pdb>

    Example:
    python create_topology.py ../5IF1/5if1_fixed_nvt1.pdb

    """
    # load in input PDB file and force field XML files
    pdb = app.PDBFile(filepath_pdf)
    force_field = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # prepare system and integrator
    system = force_field.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    integrator = mm.LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
    )
    integrator.setConstraintTolerance(0.00001)

    # Add something similar to your code.
    # It is not ideal but there seems to be an issue with parmed. It should do the job for us.
    # Add also a reference to the following github issue https://github.com/ParmEd/ParmEd/issues/930
    struct = parmed.openmm.topsystem.load_topology(pdb.topology, system, pdb.positions)
    bond_type = parmed.BondType(1.0, 1.0, list=struct.bond_types)
    struct.bond_types.append(bond_type)
    for bond in struct.bonds:
        if bond.type is None:
            bond.type = bond_type
    output_path_prefix = os.path.splitext(filepath_pdf)[0]
    struct.save(output_path_prefix + ".prmtop", overwrite=True)
    struct.save(output_path_prefix + ".crd", format="rst7", overwrite=True)
    return


def cluster_trajectory(
        filepath_topology: PathLike,
        filepath_trajectory: PathLike,
        path_output: PathLike,
        n_clusters: int = 5,
):
    """
    uses pytraj to cluster the given trajectory into a given number of clusters.

    writes summary, info, and pdb files.

    Parameters
    ----------
    filepath_topology : PathLike
        Topology file (PRMTOP format)
    filepath_trajectory
        Trajectory file (DCD) format
    path_output

    Returns
    -------

    """
    out_dir = Path(path_output)
    summary_path = str(out_dir.joinpath('cluster_summary.out'))
    info_path = str(out_dir.joinpath('cluster_info.out'))
    repout_path = str(out_dir.joinpath('repout.crd'))
    cluster_path = str(out_dir.joinpath('clusterout.crd'))

    traj = pt.load(filepath_trajectory, filepath_topology)
    frame_indices = kmeans(
        traj,
        n_clusters=n_clusters,
        mask='@CA',
        metric='rms',
        options=(
            f'summary {summary_path} info {info_path} '
            f'repout {repout_path} repfmt pdb clusterout {cluster_path}'
        )
    )
    return


def fix_pdb(filepath_pdb, chains_to_remove: List[int] = None):
    """
    Remove unwanted chains and adds missing residues in the middle of the chain.
    The result PDB file is written into a file located in the same directory as the input,
    with _fixed as suffix.

    Note that you need to pass chain indices (0, 1, ...) instead of chain IDs (A, B, ...).

    Reference:
    Script is adapted from:
    http://openmm.org/tutorials/hkmt_tip4pew/
    """
    fixer = PDBFixer(filepath_pdb)
    # remove chains if specified
    if chains_to_remove is not None:
        fixer.removeChains(chains_to_remove)
    # find missing residues
    fixer.findMissingResidues()
    # add missing residues
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    missing_residues = dict()
    for key, _ in fixer.missingResidues.items():
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missing_residues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missing_residues
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    output_path = str(Path(filepath_pdb).with_suffix("")) + "_fixed.pdb"
    with open(output_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return


prot = app.PDBFile("/Users/home/Downloads/3w32 (1).pdb")
top = prot.getTopology()