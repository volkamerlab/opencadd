# Structural Hierarchy
Many biological molecules, specifically primary metabolites, such as proteins, nucleic acids (DNA/RNA), 
and carbohydrates, are polymers and thus have a hierarchical structure: 
From a bottom-up perspective, atoms are connected together via chemical bonds, 
forming a finite set of monomers (aka residues), i.e. chemically distinct molecular fragments,
which are the building blocks of large polymeric structures; they are again connected together
via chemical bonds at defined atomic positions to form either linear or branched polymeric chains.
For example, amino acids and nucleotides, which are the monomeric building blocks of proteins and nucleic acids,
respectively, connect together to form single linear chains, while monosaccharides, the building blocks of
carbohydrates, can form both linear and branched polymers. In addition, each polymer usually interacts with
other polymers and/or small molecules via non-covalent intermolecular interactions, 
forming complex macromolecular assemblies that are responsible for various biological functions.  

The 3D structural data of the biological molecules held by the Protein Data Bank are organized 
and represented in a similar 
[hierarchical structure](https://www.rcsb.org/docs/general-help/organization-of-3d-structures-in-the-protein-data-bank) 
({ref}`Figure 1 <figure_PDB_3d_hierarchy>`), at four levels: **Entry**, **Entity**, **Instance** and **Assembly**. 
In short, an *entry* corresponds to a specific deposited structure; it contains one or more *instances* of
at least one *entity*, associated with each other in one or more *assemblies*. 
Understanding this hierarchy is essential for exploring the PDB, 
identifying relevant structures, and analyzing the acquired data.

```{figure} https://cdn.rcsb.org/rcsb-pdb/content/5fc537b03fb4b83beba83b23/Organization.jpg
:name: figure_PDB_3d_hierarchy

Figure 1. Hierarchical organization of 3D structures in PDB. 
**A)** The *entry* 2hbs contains two complete sickle cell hemoglobin tetramers 
(the two large blobs on the top-right and bottom-left corners), each of which in complex with several 
heme cofactors (red blobs), and surrounded by water molecules (small green spheres).
**B)** Each tetramer (one in grey, and the other in original colors) is an *assembly*, 
i.e. the functional biological unit in cells, responsible for binding to oxygen and its delivery in the blood.
Despite being nearly identical, each *assembly* in the *entry* has a unique Biological Assembly ID in the PDB databse. 
Moreover, each assembly consists of two polymer *entities* (aka chains): alpha and beta. 
More specifically, each assembly contain four *instances* of polymer *entities*: two *instances* of 
the alpha *entity* (in organe and yellow), and two *instances* of the beta *entity (in two shades of blue)*. 
In addition, each assembly includes four *instances* of the heme cofactor 
(each associated with one of the four polymer *instances*), and several *instances* of water.
**C)** The heme cofactor (in red) and water (in green) are non-polymeric *entities*. The whole *entry* contains 
eight *instances* of heme (each bound to one *instance* of one polymer *entity*), plus several hundred *instances*
of water.
**D)** All *instances* of all polymeric *entities* are colored, while heme and water entities are shown in grey.
**E)** All four *instances* of the alpha *entity* are colored, while all other *entities* are shown in grey
(*instances* of the water *entity* are not shown for clarity).
**F)** Only one *instance* of the alpha *entity* is highlighted.
Image retrieved from the [RCSB PDB website](https://www.rcsb.org/docs/general-help/organization-of-3d-structures-in-the-protein-data-bank).
```

## Entry
An *entry* corresponds to all the data for a particular 3D structure deposited in the PDB. It is identified 
by a 4-character alphanumeric code, called the PDB ID (e.g. 2hbs). 
Every *entry* in the PDB archive contains at least one *entity*.

## Entity
An *entity* is a chemically unique molecule; 
it can be polymeric (e.g. a single protein chain or nucleic acid strand), 
branched (e.g. a polysaccharide), or non-polymeric (e.g. a ligand or a solvent molecule).

```{note}
In the PDB, oligo-/polysaccharides are called *branched entities*, regardless of whether their actual 
polymeric structure is linear or branched.
```

## Instance
An *entry* may contain multiple identical *entities*, i.e. multiple occurrences of the same molecule.
Each occurrence of a particular *entity* is called an *instance* of that *entity*.
In case of polymeric and branched *entities* (aka chains), each *instance* is given a unique alphanumeric chain ID (e.g. A, A1, AC1 etc.).
On the other hand, *instances* of non-polymeric *entities* are identified by the chain ID of the closest 
neighboring *instance* of a polymeric *entity*, in addition to a unique sequence number. 
For example, two heme groups associated with the same protein chain *instance* **A** 
may be identified as **A101** and **A102**.

```{note}
Within an *entry*, chain IDs can be used to uniquely identify each *instance* of each polymeric/branched *entity*.
However, assignment of chain IDs is not systematic throughout the PDB archive, meaning that the same *entity*
in two different PDB *entries* may have different chain IDs assigned to its *instances*.
```

## Assembly
An *assembly* is a biologically relevant group of one or more *instances* of one or more *entities* that are
associated with each other to form a stable complex and/or perform a specific function. Each *assembly* is
identified by a unique Biological Assembly ID, which is a positive integer starting from 1. 
