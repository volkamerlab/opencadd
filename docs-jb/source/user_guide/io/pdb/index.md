---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# PDB [`opencadd.io.pdb`](../../../api_reference/_autosummary/opencadd.io.pdb.rst)

## Background Overview 
The Protein Data Bank (PDB) file format is a commonly used biochemical data file,
containing information on biochemical systems, such as nucleic acids, proteins, protein–ligand complexes, 
and biologically relevant small molecules. \
Learn more about the Protein Data Bank and its file format on the next page: [](theory.md)

## Usage Overview
openCADD is able to fully process almost every bit of information in a PDB file, and present in a concise 
data structure; this makes it very easy to inspect, analyze, and manipulate the data, which is usually the 
next step after importing. The parser is also faster than most other alternatives, and is even able to 
selectively parse only specific parts of the file, for more speed gains (this is specially useful when parsing
a large number of files, where only a specific set of information is of interest).\
Learn more about 


## TL;DR
Import `opencadd`:
```{code-cell} ipython3
import opencadd as oc
```
and load a PDB file, for example using its PDB ID (requires internet connection):
```{code-cell} ipython3
pdb_3w32 = oc.io.pdb.from_pdb_id("3w32")
```
or from a local PDB file on your device:
```{code-cell} ipython3
pdb_3w32 = oc.io.pdb.from_filepath("my_files/my_copy_of_3w32.pdb")
```
For more options, see [Usage](usage.md).

Every piece of information in the PDB file is now at your fingertips:




```{code-cell} ipython3
:tags: [mytag]

import opencadd as oc
pdb_3w32 = oc.io.pdb.from_pdb_id("3w32", parse=False)
pdb_3w32.records_atom_hetatm
```





```{toctree}
:maxdepth: 2
:hidden:

theory
usage
examples
```