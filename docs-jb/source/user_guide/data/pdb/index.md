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

# PDB `opencadd.db.pdb`
The [Protein Data Bank](https://www.wwpdb.org/) (PDB) is the largest database 
of experimentally determined structural data for large biological molecules, 
such as proteins, nucleic acids, and related complex assemblies. 

```{admonition} Learn More
Learn more about the Protein Data Bank in the Learn section: [](../../../learn/databases/pdb/index.md) 
```

## Usage
PDB provides two main services for programmatic access to its database:
* **File Download Services** allows for downloading various data files from the PDB archive, 
  such as PDB entry files, small molecules files, experimental data files and electron density maps etc.
* **Web Services** including a search API and a data API, to search the database for complex queries 
  and obtain specific data.

Currently, `opencadd.db.pdb` provides access to the REST data API and file download services of the PDB.

First, import `opencadd`:
```{code-cell} ipython3
import opencadd as oc
```

### File Download Services `opencadd.db.pdb.file`
PDB entry files and small molecule files can be downloaded in various formats.

#### PDB Entry Files
A PDB entry file can be downloaded by providing its PDB ID to the `entry` function.
For example, to download the PDB entry 3W32:

```{code-cell} ipython3
pdb_3w32_asym_cif = oc.db.pdb.file.entry("3w32")
```

```{note}
The PDB ID in this case is not case sensitive, i.e. "3W32" and "3w32" are equivalent.
```

By default, this downloads the asymmetric unit of the PDB entry. If instead a certain biological assembly
of the PDB entry is needed, it can be downloaded by providing its Biological Assembly ID.

```{admonition} Learn More
Learn more about the differences between asymmetric units and biological assemblies
in the Learn section: [](../../../learn/databases/pdb/index.md)
```

For example, to download the assembly 1 of the PDB entry 3W32:

```{code-cell} ipython3
pdb_3w32_assembly1_cif = oc.db.pdb.file.entry(
  "3w32", 
  biological_assembly_id=1
)
```

```{tip}
Biologicall Assembly IDs for each entry start at 1 and increase by 1.
```

By default, all entry files are downloaded in the PDBx/mmCIF file format. To download the entry file in 
another format, the required format must be specified. 
For example, to download the same assembly in PDB format:

```{code-cell} ipython3
pdb_3w32_assembly1_pdb = oc.db.pdb.file.entry(
  "3w32", 
  biological_assembly_id=1, 
  file_format="pdb"
)
```

In general, four different file formats can be chosen: "cif" (default), "pdb", "xml", and "bcif".
However, it should be noted that some PDB entries are not provided in the PDB format 
([learn more](../../../learn/databases/pdb/index.md)).

By default, the function returns the content of the downloaded entry file in `bytes`:

```{code-cell} ipython3
type(pdb_3w32_assembly1_pdb)
```

Let's decode the bytes content to `str` and take a look at the first 400 characters:
```{code-cell} ipython3
print(pdb_3w32_assembly1_pdb.decode()[:400])
```

Another option is to specify a local path in order to directly save the file on your device:
```{code-cell} ipython3
pdb_3w32_assembly1_pdb_filepath = oc.db.pdb.file.entry(
  "3w32", 
  biological_assembly_id=1, 
  file_format="pdb",
  output_path="./my_pdb_archive"
)
```
In this case, the function returns the full path of the stored file as a `pathlib.Path` instance:
```{code-cell} ipython3
type(pdb_3w32_assembly1_pdb_filepath)
```
The filename will always be the PDB ID of the downloaded entry:
```{code-cell} ipython3
pdb_3w32_assembly1_pdb_filepath.name
```

```{hint}
openCADD has a powerful PDB parser, which can process all data in a PDB file and present it
in a useful data structure for further analysis. This is discussed in the next chapter: 
[](../../io/pdb/index)
```

#### Small Molecule Files
In addition to PDB entry files, small molecule files can also be downloaded in various file formats 
using the `small_molecule` function. Small molecules include ligands and chemical components
maintained in the Chemical Component Dictionary (CCD) and the Biologically Interesting Molecule 
Reference Dictionary (BIRD). These are identified, and can be downloaded, by their Ligand ID.
For example, to download the ligand W32 (the co-crystallized ligand in PDB entry 3W32):

```{code-cell} ipython3
lig_w32_model_sdf = oc.db.pdb.file.small_molecule("w32")
```

For each ligand, three types of data files can be obtained: model coordinates, ideal coordinates,
and definition data. By default, model coordinates are downloaded. In order to download ideal coordinates or
definition files, 'ideal_coord' and 'def' can be provided as the argument of the `file_type` parameter, respectively.
For example, to download ideal coordinates of the same ligand:

```{code-cell} ipython3
lig_w32_ideal_sdf = oc.db.pdb.file.small_molecule(
  "w32",
  file_type="ideal_coords",
)
```

While definition files are only available in CIF format, coordinates data can be downloaded 
in either SDF or MOL2 file formats. For these files, the SDF format is downloaded by default.
In order to download a coordinates file, e.g. the ideal coordinates of the above ligand, in MOL2 format, 
argument 'mol2' must be provided to the `file_format` parameter:

```{code-cell} ipython3
lig_w32_ideal_mol2 = oc.db.pdb.file.small_molecule(
  "w32",
  file_type="ideal_coords",
  file_format="mol2",
)
```

Similar to PDB entry files, the content is returned in `bytes` by default, but a local path can be provided
to the `output_path` parameter to instead directly store the file on your local device and get the 
full filepath as the return value:

```{code-cell} ipython3
lig_w32_ideal_mol2_filepath = oc.db.pdb.file.small_molecule(
  "w32",
  file_type="ideal_coords",
  file_format="mol2",
  output_path="./my_lig_archive"
)
```

Again, similar to PDB entry files, the filename will always be the Ligand ID:
```{code-cell} ipython3
lig_w32_ideal_mol2_filepath.name
```


### Data API `opencadd.db.pdb.data`
Various data at different levels can be obtained.

#### Repository Holdings Data
Information about entry holdings in the PDB archive can be retrieved, using functions that start with `holdings`.
For example, to get all PDB IDs of current entries in the PDB archive:
```{code-cell} ipython3
all_curr_pdb_ids = oc.db.pdb.data.holdings()
```
This returns a list of PDB IDs as a `numpy` array:
```{code-cell} ipython3
all_curr_pdb_ids
```
The total number of current PDB entries in the PDB archive at this moment
is thus equal to the size of this array:
```{code-cell} ipython3
import datetime

print(
  f"As of {datetime.date.today()}, the PDB archive contains "
  f"{all_curr_pdb_ids.size} entries."
)
```
Instead of all current entries, the PDB IDs of all "unreleased" and "removed" entries can also be obtained:
```{code-cell} ipython3
all_unreleased_pdb_ids = oc.db.pdb.data.holdings("unreleased")
all_removed_pdb_ids = oc.db.pdb.data.holdings("removed")
```

As mentioned above, some PDB entries are only available in mmCIF file format, and do not have corresponding
PDB files; the PDB IDs of all such entries can be obtained:
```{code-cell} ipython3
oc.db.pdb.data.holdings_without_pdb_file()
```

Inversely, given a PDB ID, its holdings status can be checked using the `holdings_status` function:

```{code-cell} ipython3
oc.db.pdb.data.holdings_status("3w32")
```

```{attention}
Currently, this and most other functions in `opencadd.db.pdb.data` return data in form of nested dictionaries,
in compliance with the PDBx/mmCIF Dictionary schema. This may change to a more efficient data structure 
in the future.
```


In addition, description of a specific unreleased or removed entry can be retrieved, providing its PDB ID to
the `holdings_unreleased` and `holdings_removed` functions, respectively:
```{code-cell} ipython3
some_unreleased_pdb_id = all_unreleased_pdb_ids[0]
oc.db.pdb.data.holdings_unreleased(some_unreleased_pdb_id)
```
```{code-cell} ipython3
some_removed_pdb_id = all_removed_pdb_ids[0]
oc.db.pdb.data.holdings_removed(some_removed_pdb_id)
```

#### Groups Level Data

#### Entry Level Data

#### Assembly Level Data

#### Entity Level Data

#### Instance Level Data

#### Chemical Component Data

#### Interface Data

#### Dictionary Schema

```{toctree}
:maxdepth: 2
:hidden:

```

