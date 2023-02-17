
`wwPDB Documentation <https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html>`_

## Section: Title
contains a wide range of information about the entry as a whole,
including a title describing the experiment, identifiers, classification, experimental techniques, 
biological molecules, corresponding publications and authors.
This section contains the largest number of records among all sections; these are:
HEADER, OBSLTE, TITLE, SPLIT, CAVEAT, COMPND, SOURCE, KEYWDS,
EXPDTA, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and JRNL.

### Record: HEADER
is a mandatory record that must appear as the first record
of every PDB file. It contains three pieces of information, corresponding to three fields:

#### Field: idCode
contains the PDB identification code (PDB ID) of the entry.
This is a 4-character unique identifier within the Protein Data Bank (PDB), e.g. 3W32.
The first character must be a digit in the range 0 - 9, 
while the remaining three characters are alphanumeric, i.e. may either be digits or letters.
```{info}
Entries whose PDB ID starts with 0 do not contain coordinate data;
these so-called 'no coordinates' (NOC) files are now removed from the PDB archive.
```
#### Field: classification
contains a set of classification tags for each biologically interesting molecule 
present in the entry. The classifiers may be based on, e.g. function, metabolic role, 
molecule type, cellular location etc., but according to 
[Section B](https://cdn.rcsb.org/wwpdb/docs/documentation/annotation/wwPDB-B-2022Nov-V4.2-2.pdf)
of the [wwPDB Deposition and Biocuration Documentation](https://www.wwpdb.org/documentation/biocuration),
they must match a set of pre-defined categories. This means that each tag must either exactly match one of the
classification tags [listed here](https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/class.dat), 
or it can be a compound word made by adding words such as INHIBITOR, ACTIVATOR,
RECEPTOR, REGULATOR etc. to one of the listed classification tags.
For example, the TRANSFERASE tag can be expanded to TRANSFERASE INHIBITOR, TRANSFERASE ACTIVATOR etc.
```{note}
Due to the limited length of the classification field, tags must sometimes be abbreviated. 
In such cases, the full terms also appear in KEYWDS records, in no strict order.
```
#### Field: depDate
contains the date of deposition of the entry at the Protein Data Bank.

### Record: OBSLTE
appears only in entries that have been “obsoleted”, i.e. removed from public distribution, 
and possibly replaced by a more accurate entry.
This is done when the structure's coordinates undergo major revisions, 
resulting in a change of the structure's geometry or chemical composition, 
such as changes in polymer sequences, or identity of ligands.
```{info}
While having been removed from the main PDB distribution, all OBSLTE entries are still 
separately accessible via the [PDB archive](https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete>).
```
Similar to the HEADER record, this record also contains three pieces of data in three fields:
#### Field: idCode
contains the PDB ID of the entry; same as in the HEADER record.
#### Field: repDate
contains the date of replacement/removal of the entry from the PDB distribution.
#### Field: rIdCode
contains a list of PDB IDs of the new entries that have replaced (i.e. superseded) this entry.
A reference to the obsoleted entry can be found in all superseding entries.

### Record: TITLE
is a mandatory record containing a free-form text, describing the experiment, analysis, 
and/or contents that is represented in the entry, 
and any procedures or conditions that distinguish it from similar entries.
Some data that may be included are experiment type, description of the mutation,
and the fact that only alpha carbon coordinates have been provided in the entry.

### Record: SPLIT
is only present in entries that constitute a part of a larger entry.
As discussed earlier, the PDB file format can only contain coordinates for 99,999 atoms, 
and no more than 62 chains. Therefore, larger entries that exceed these limitations 
must be split into multiple PDB files. In such cases, the SPLIT record contains
a list of all PDB IDs of entries that are required to reconstitute the complete complex.
```{hint}
A related piece of information can also be found in REMARK 350 records,
describing the entire complex.
```

### Record: CAVEAT
appears only in entries containing errors or unresolved issues, 
or when the wwPDB is unable to verify the transformation of the coordinates 
back to the crystallographic cell. Notice that in the latter case, 
the molecular structure may still be correct. CAVEAT has two fields:
#### Field: idCode
repeats the PDB ID of the entry; same as in HEADER record.
#### Field: comment
contains a free-form text description of the caveat.

### Record: COMPND
is a mandatory record that describes either the macromolecules in the entry,
or a standalone drug or inhibitor in cases where the entry does not contain a polymer.
COMPND contains a single field, also named *compound*.
#### Field: compound
contains a list of specifications for each (macro)molecule in the entry; these are:

* MOL_ID: Enumerates each molecule; the same ID appears also in the SOURCE records.
* MOLECULE: Name of the (macro)molecule. 
  For chimeric proteins, the protein name is 
  comma-separated and may refer to the presence of a linker, e.g. "protein_1, linker, protein_2".
* CHAIN: List of chain IDs labeling the instances of this molecule present in the entry. 
* FRAGMENT: Name of a specific domain or region in the molecule. 
* SYNONYM: List of synonyms for the name of the molecule, i.e. synonyms for MOLECULE.
* EC: List of Enzyme Commission numbers for the molecule.
* ENGINEERED: Whether the molecule is engineered, i.e. produced using recombinant technology or chemical synthesis.
* MUTATION: Whether the macromolecule contains a mutation.
* OTHER_DETAILS: Additional free-form text comments.

```{note}
Other than MOL_ID, MOLECULE and CHAIN, it is not mandatory for other specifications 
to be listed for every molecule present in the entry. Moreover, one molecule (with the same MOL_ID and
MOLECULE) may appear multiple times, each time corresponding to a ceratin FRAGMENT within the molecule.
```

### Record: SOURCE
is a mandatory record that contains information on the biological/chemical source of
each biological molecule in the PDB file, or a standalone drug or inhibitor in cases
where the entry does not contain a polymer. Like the COMPND record, SOURCE also contains a 
single field: *srcName*
#### Field: srcName
is similar to the *compound* field of the COMPND record, and contains a specification list for
each macromolecule in the entry, with following information:
* MOL_ID: Enumerates each molecule; the same ID appears also in the COMPND records.
* SYNTHETIC: Whether the molecule was prepared purely by chemical synthesis; 
  this may be indicated by a simple 'YES', or by phrases such as "NON-BIOLOGICAL SOURCE" or 
  "BASED ON THE NATURAL SEQUENCE". The `engineered` column in the COMPND record 
  is also set in such cases.
* FRAGMENT: Name of a specific domain or region in the molecule, like in COMPND records. 
* ORGANISM_SCIENTIFIC: Scientific name of the organism.
* ORGANISM_COMMON: Common name of the organism.
* ORGANISM_TAXID: NCBI Taxonomy ID of the organism.
* STRAIN: Identifies the strain.
* VARIANT: Identifies the variant.
* CELL_LINE: The specific line of cells used in the experiment.
* ATCC: American Type Culture Collection tissue culture number.
* ORGAN: Organized group of tissues that carries on a specialized function.
* TISSUE: Organized group of cells with a common function and structure.
* CELL: Identifies the particular cell type.
* ORGANELLE: Organized structure within a cell.
* SECRETION: Identifies the secretion, such as saliva, urine, or venom, from which the molecule was isolated.
* CELLULAR_LOCATION: Identifies the location inside/outside the cell, where the compound was found. Examples are: 'extracellular', 'periplasmic', 'cytosol'.
* PLASMID: Identifies the plasmid containing the gene.
* GENE: Identifies the gene.
* EXPRESSION_SYSTEM: Scientific name of the expression system; i.e. name of the organism in which the molecule was expressed.
* EXPRESSION_SYSTEM_COMMON: Common name of the expression system. 
* EXPRESSION_SYSTEM_TAXID: NCBI Taxonomy ID of the expression system.
* EXPRESSION_SYSTEM_STRAIN: Strain of the organism in which the molecule was expressed.
* EXPRESSION_SYSTEM_VARIANT: Variant of the organism used as the expression system.
* EXPRESSION_SYSTEM_CELL_LINE: The specific line of cells used as the expression system.
* EXPRESSION_SYSTEM_ATCC_NUMBER: American Type Culture Collection tissue culture number of the expression system.
* EXPRESSION_SYSTEM_ORGAN: Specific organ which expressed the molecule.
* EXPRESSION_SYSTEM_TISSUE: Specific tissue which expressed the molecule.
* EXPRESSION_SYSTEM_CELL: Specific cell type which expressed the molecule.
* EXPRESSION_SYSTEM_ORGANELLE: Specific organelle which expressed the molecule.
* EXPRESSION_SYSTEM_CELLULAR_LOCATION: Identifies the location inside or outside the cell which expressed the molecule.
* EXPRESSION_SYSTEM_VECTOR_TYPE: Identifies the type of vector used, i.e. plasmid, virus, or cosmid.
* EXPRESSION_SYSTEM_VECTOR: Identifies the vector used.
* EXPRESSION_SYSTEM_PLASMID: Plasmid used in the recombinant experiment.
* EXPRESSION_SYSTEM_GENE: Name of the gene used in recombinant experiment.
* OTHER_DETAILS: Other details about the source as free-form text.

```{note}
Similar to COMPND records, only the relevant specifications for each molecule are provided.
Moreover, when necessary to fully describe hybrid molecules, they may be listed multiple times.
Hybrid molecules prepared by fusion of genes are treated as multi-molecular systems for
the purpose of specifying the source. Similar to COMPND records, FRAGMENT is used to associate the data
with its corresponding fragment.
```

### Record: KEYWDS
is a mandatory record, containing a list of free-form keywords/terms relevant to the PDB file, 
similar to that found in journal articles. The keywords may for example describe 
functional classification, metabolic role, known biological or chemical activity, 
or structural classification.
```{note}
All classification tags present in the 
*classification* field of the HEADER record are also repeated here, unabbreviated. 
However, unlike the *classification* field where the tags are grouped per molecule, 
here all keywords are simply listed in no order.
```

### Record: EXPDTA
is a mandatory record, containing a list of pre-defined tags to identify
the experimental techniques used for determining the structure. The allowed tags are:

* X-RAY DIFFRACTION
* FIBER DIFFRACTION
* NEUTRON DIFFRACTION
* ELECTRON CRYSTALLOGRAPHY
* ELECTRON MICROSCOPY
* SOLID-STATE NMR
* SOLUTION NMR
* SOLUTION SCATTERING

### Record: NUMMDL

### Record: MDLTYP

### Record: AUTHOR

### Record: REVDAT

### Record: SPRSDE

### Record: JRNL

## Section: Remarks

## Section: Primary Structure 
Primary structure section of the PDB file, with information on the primary structure of the polymers
present in the PDB file, including the sequence of residues in each chain of the macromolecule(s);
references to corresponding sequences in other databases, and possible conflicts between the two;
and modifications of the residues in each sequence.

* This section contains all data from the PDB records DBREF, DBREF1, DBREF2, SEQADV, SEQRES, and MODRES.
* Embedded in these records are chain identifiers and sequence numbers that allow other records
  to link into the sequence.


## Section: Heterogen
Heterogen section of the PDB file, containing the complete description of non-standard
residues in the entry, such as prosthetic groups, inhibitors, solvent molecules, and ions,
for which coordinates are supplied.

Notes
-----
* Groups are considered HET if they are not part of a biological polymer described in
  SEQRES and considered to be a molecule bound to the polymer, or they are a chemical
  species that constitute part of a biological polymer and is not one of the following:
  * standard amino acids
  * standard nucleic acids (C, G, A, U, I, DC, DG, DA, DU, DT and DI)
  * unknown amino acid (UNK) or nucleic acid (N) where UNK and N are used to indicate
    the unknown residue name.
* HET records also describe chemical components for which the chemical identity is
  unknown, in which case the group is assigned the hetID UNL (Unknown Ligand).
* A particular HET group is represented in the PDB archive with a unique hetID.
* PDB entries do not have HET records for water molecules, deuterated water,
  or methanol (when used as solvent).
* Unknown atoms or ions will be represented as UNX with the chemical formula X1.
  Unknown ligands are UNL; unknown amino acids are UNK.

References
----------
* Detailed chemical definitions of non-polymer chemical components are described in the
  Chemical Component Dictionary at: https://ftp.wwpdb.org/pub/pdb/data/monomers


## Secondary Structure Section
Secondary structure section of the PDB file, describing helices and sheets found in the
protein and polypeptide structures.


## Miscellaneous Features Section
Miscellaneous features section of the PDB file, describing properties in the molecule such as
environments surrounding a non-standard residue or the assembly of an active site.


## Connectivity
Connectivity and Connectivity Annotation sections of the PDB file, providing information on atomic connectivity.


## Bookkeeping Section
Bookkeeping section of a PDB file. It contains the information from the MASTER record,
which is a control record, counting the occurrence (checksum) of selected record types
in the PDB file.