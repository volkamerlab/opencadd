
# Finding Data [`opencadd.db`](../../api_reference/_autosummary/opencadd.db.rst)
Data plays a critical role in many CADD and bioinformatics studies. 
Most pipelines rely heavily on some form of input data, 
such as structural data of biologically active molecules; 
sequence data of biological polymers, e.g. nucleic acids and proteins; 
or interaction data, e.g. the half maximal inhibitory concentration (IC{sub}`50`). 
A typical first step in a CADD pipeline is therefore finding and gathering all the data required
for the specific project at hand.

```{admonition} Learn More
Learn more about CADD-relevant data and databases in the Learn section: [](../../learn/databases/index.md) 
```

The aim of the `opencadd.db` package is to provide full programmatic access 
to a variety of open access online databases containing biochemical and pharmaceutical data. 
It allows users to query different databases, 
to find and retrieve the exact set of data they need for their project. 
These can be in the form of various file formats and data structures, which can either be directly worked with, 
or used as input in the next steps of the pipeline.

Currently, following databases can be accessed:

* Protein Data Bank (PDB) `opencadd.db.pdb`
* Kinase–Ligand Interaction Fingerprint and Structure (KLIFS) database `opencadd.db.klifs`
* PubChem `opencadd.db.pubchem`


```{toctree}
:maxdepth: 2
:hidden:

pdb/index
pubchem/index
embl/index
klifs/index
```
