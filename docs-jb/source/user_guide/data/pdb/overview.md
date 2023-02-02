# Overview

## History of Protein Data Bank
The development of experimental techniques such as X-ray crystallography and nuclear magnetic resonance (NMR)
led to substantial advancements in the field of structural biology. A great breakthrough was the successful 
modeling of the DNA double-helix structure, by Watson and Crick in 1953, based on previous research by
Rosalind Franklin and Maurice Wilkins.{cite}`Watson&Crick` Later in 1958, John Kendrew determined the 
tertiary structure of a protein, myoglobin, for the first time.{cite}`first_determined_tertiary_protein_structure` 
These discoveries catapulted the field of structural biology to a new era, and led to a revolution in 
numerous related fields, such as medicinal chemistry, pharmacology, and medicine. 

In the days before the internet, researchers around the world needed a way to share the data on these
experimentally determined macromolecular crystal structures, which were quiet large, 
containing the three-dimensional coordinates of thousands of atoms for each macromolecule.
To address this issue, the Protein Data Bank (PDB) was established in 1971 at 
Brookhaven National Laboratories (BNL), as the first open access digital data resource 
in biology and medicine.{cite}`PDB_first_announcement, PDB_first_release` 
Structural biologists would send their data to the PDB, who would then store them on magnetic tapes, 
and send them to interested users via postal mail, for no charge other than the mailing fees. 
Having started with only a handful of available crystal structures ({ref}`Figure 1 <figure_first_pdb_structures>`), 
the size of the archive began to grow rapidly in the early 1990s ({ref}`Figure 2 <pdb_growth>`), 
thanks to the improvement of experimental techniques, among others.{cite}`PDB_history_review`

```{figure} https://cdn.rcsb.org/rcsb-pdb/v2/about-us/early.png
:name: figure_first_pdb_structures

Figure 1. The full archive of Protein Data Bank (PDB) in 1973 consisted of only nine macromolecular structures. 
Image retrieved from the [RCSB PDB website](https://www.rcsb.org/pages/about-us/history).
```

In 1998, the Research Collaboratory for Structural Bioinformatics (RCSB) became responsible for the 
management of the PDB.{cite}`RCSB_PDB_first_release` Shortly afterwards, 
the Macromolecular Structure Database (MSD) at the European Bioinformatics Institute (EBI) 
and the Protein Data Bank Japan (PDBj) at the Institute for Protein Research in Osaka University
joined the U.S. based RCSB to form the worldwide Protein Data Bank (wwPDB) in 2003, with the goal of
maintaining a single freely and publicly available archive of macromolecular structural data for 
the global community.{cite}`wwPDB_first_release`

```{figure} images/growth_of_PDB.jpg
:name: pdb_growth

Figure 2. Evolution of the PDB archive; aggregated (A) and individual (B) number of released structures per year, 
resolved and colored by experimental technique (MX: Macromolecular X-ray Crystallography, NMR: Nuclear Magnetic 
Resonance, 3DEM: 3D Electron Microscopy).{cite}`PDB_report2018` Image retrieved from reference {cite}`PDB_report2018`.
```

Today, the PDB is the largest and most comprehensive database of its kind, hosting detailed experimental data 
on more than 214,000 (as of February 2023) structures. The cost to replicate the contents of the PDB archive 
is now estimated at $18 billion U.S. Dollars. (more statistics available on 
[RCSB](https://www.rcsb.org/stats/) and [wwPDB](http://www.wwpdb.org/stats/deposition)).
Since understanding the structures of biological macromolecules and their interactions is a crucial step in
many research fields such as drug development, the PDB has become one of the most essential resources for 
both computational and experimental researchers working on related subjects. Almost all major scientific 
journals now require researchers to submit their structure data to the PDB for their publications.

The full PDB archive is freely accessible, and can be programmatically accessed via a provided web API. 

# References
```{bibliography}
:style: unsrt
```