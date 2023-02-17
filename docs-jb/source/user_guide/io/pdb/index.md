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

The Protein Data Bank (PDB) file format is a commonly used biochemical data file,
containing information on biochemical systems, such as nucleic acids, proteins, protein–ligand complexes, 
and biologically relevant small molecules.

```{admonition} Learn More
Learn more about the PDB file format in the Learn section: [](../../../learn/databases/pdb/file_format.md) 
```

## Usage
openCADD is able to fully process almost every bit of information in a PDB file, and present in a concise 
data structure, making it very easy to inspect, analyze, and manipulate the data.

First, import `opencadd`:
```{code-cell} ipython3
import opencadd as oc
```

### Importing PDB Files
PDB files can be imported and parsed from different sources:

#### From a Local File
PDB files stored locally on your device can be imported using the `from_filepath` function:
```{code-cell} ipython3
my_pdb_filepath = oc.db.pdb.file.entry(
  "3W32", 
  file_format="pdb",
  output_path="./my_pdb_files"
)
my_pdb = oc.io.pdb.from_filepath(my_pdb_filepath)
```
```{note}
In the above example, the `opencadd.db` package is used to download and store the PDB file
of the entry 3W32 in a local directory named 'my_pdb_files'. The full filepath of the downloaded file
is stored in the `my_pdb_filepath` variable, which is then used to import the PDB file using the 
`from_filepath` function. In a real scenario, the user already has a PDB file stored locally on their device,
and can directly use the corresponding filepath to import the file using the `from_filepath` function.
```

#### From File Content
It may be that the content of the PDB file has already been retrieved 
in some way during the runtime; in such cases, instead of having to first write the data to a file and then
use the `from_filepath` function to load the PDB file, the file content in memory can be directly
inputted to the `from_file_content` function:
```{code-cell} ipython3
my_pdb_content = oc.db.pdb.file.entry(
  "3W32", 
  file_format="pdb",
)
my_pdb = oc.io.pdb.from_file_content(my_pdb_content)
```
```{note}
Again, we are using the `opencadd.db.pdb.file.entry` function to obtain the PDB file content, 
for demonstration purposes; notice that in contrast to the first example, this time no argument is
provided to the `output_path` parameter of the function. This will result in the function returning the
downloaded PDB file's content directly, instead of storing it on disk and returning the filepath.
However, in a real scenario, the user already has the PDB file's content in memory,
and can directly use the `from_file_content` function.
```

#### From PDB ID
When the desired PDB file is a standard PDB entry available in the online PDB archive,
the most convenient way is to provide the entry's PDB ID to the `from_pdb_id` function:
```{code-cell} ipython3
pdb_3w32 = oc.io.pdb.from_pdb_id("3w32")
```
By default, this will download and parse the asymmetric unit of the PDB entry. 
If instead a certain biological assembly of the PDB entry is needed, 
it can be downloaded by providing its Biological Assembly ID.

```{admonition} Learn More
Learn more about the differences between asymmetric units and biological assemblies
in the Learn section: [](../../../learn/databases/pdb/index.md)
```

```{attention}
Some data are only available in the asymmetric unit files; if these data are needed in addition to
the coordinates of a specific biological assembly, you can import both files and pull the relevant
data from each version. 
```

For example, to fetch and parse the assembly 1 of the PDB entry 3W32:

```{code-cell} ipython3
pdb_3w32_assembly1 = oc.io.pdb.from_pdb_id(
  "3w32",
  biological_assembly_id=1,
)
```

```{tip}
Biologicall Assembly IDs for each entry start at 1 and increase by 1.
```

```{hint}
This is a convinience function that uses the `opencadd.db.pdb.file.entry` function 
to download the PDB file and load it onto memory, and subsequently calls the `from_file_content` function,
just as demonstrated in the [](#### From File Content) section above.
```

### Accessing Data
Regardless of the way the PDB file was imported, all above functions fully process the data and return it as
an `opencadd.io.pdb.datastruct.PDBFile` instance. This is a data structure that closely resembles the
structure of a PDB file, allowing for direct access to all data.

```{admonition} Learn More
Learn more about the structure of PDB files and the different types of information they contain 
in the Learn section: [](../../../learn/databases/pdb/file_format.md) 
```


### Selective Parsing
The parser is also faster than most other alternatives, and is even able to 
selectively parse only specific parts of the file, for more speed gains (this is specially useful when parsing
a large number of files, where only a specific set of information is of interest).\
Learn more about 

### Data Validation
By default, the parser assumes that the input PDB 



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