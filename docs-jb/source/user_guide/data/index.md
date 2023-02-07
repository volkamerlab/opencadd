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

# Finding Data [`opencadd.data`](../../api_reference/_autosummary/opencadd.data.rst)

## TL;DR
`opencadd.data` allows for full programmatic access to several online databases, including the 
Protein Data Bank (PDB), European Molecular Biology Laboratory (EMBL), PubChem, and the 
Kinase–Ligand Interaction Fingerprint and Structure (KLIFS) database.

The following are a few examples of how `opencadd.data` works. 

First, import `opencadd`:
```{code-cell} ipython3
import opencadd as oc
```


The [`opencadd.data`](../../api_reference/_autosummary/opencadd.data.rst) package allows users to read/write.



Currently, the following formats are supported:

## Chemical Structure Data
Chemical structure data


```{toctree}
:maxdepth: 2
:hidden:

pdb/index
pubchem/index
embl/index
klifs/index
```