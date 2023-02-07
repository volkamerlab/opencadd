# PDB File Format

The PDB file format was invented by the Protein Data Bank in 1976, as a textual file format
for sharing data on experimentally determined macromolecular structures.{cite}`PDB_history_review`
It was designed for human-readability, and compatibility with [punch cards](https://en.wikipedia.org/wiki/Punched_card),
which were the common medium for digital data storage in those days. As a result, many compromises had
to be made with regard to machine-readability. After numerous revisions, version 3.30, initially released on July 2011, 
was announced to be the last version of the PDB file format, with no future changes being made. 
Later in 2014, the Macromolecular Crystallographic Information (mmCIF) 
file format replaced PDB as the standard format for the PDB archive. In 2019, 
the wwPDB [announced](https://www.wwpdb.org/news/news?year=2019#5d17c1ceea7d0653b99c87e2) that crystallographic
depositions would only be accepted in the mmCIF file format. 

Nevertheless, the wwPDB continues to accept NMR and EM data in PDB format, 
and provides PDB files (compliant with version 3.30) for all entries 
(regardless of the deposition format), unless not possible due to the limitations of the PDB file format, 
for example in case of large macromolecular complexes. In addition, PDB remains to be a common file format 
in the bioinformatics ecosystem, and many software packages still use it as their primary input/output format. 

## PDB File Format
https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf

### Entries Without PDB Format Files
The PDB file format is unable to hold data for large macromolecular structures, 
specifically, structures with more than 99,999 atoms, or with more than 62 chains.
This is due to the fact that the PDB format has a fixed column size and data type for tabular contents;
The column for atom serial numbers has a width of 5 characters, and thus can only hold integer values
between 1 and 99,999. Similarly, the chain ID column is assigned only a single alphanumerical character, 
and can only have 62 unique values (i.e. 26 lower-case leters, 26 upper-case letters, and 10 digits).

Before the phase-out of the PDB file format in 2014, these structures were 'split' into several PDB files,
and the complete structure had to be re-assembled by the user from a collection of split files. After the
replacement of the PDB file format by the PDBx/mmCIF format, these files were combined into single mmCIF files, 
and the corresponding PDB files are now only available only in some cases, containing only coordinates data.
For more details, refer to [RCSB Documentation: Structures Without Legacy PDB Format Files](https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files).


There are also some other added structural data that do not fit the old PDB file format. 
In summary, the Protein Data bank does not offer legacy PDB-format files for:
* entries with multiple-character chain IDs
* entries with more than 62 chains
* entries with 100,000 or more atoms
* entries with a complex beta sheet topology




### Contents of PDB Files

The main piece of information in a PDB file is the three-dimensional structural information on 
the atoms in the corresponding system (e.g. a nucleic acid, or a protein-ligand complex in water), 
including their positions, element types, and the residues/molecules they belong to. 
In addition, PDB files contain many other types of information, including details on the molecules in the system, 
such as classification, sources, primary and secondary structure information; 
details on the experiment and the techniques used;
database references, and even bibliographic data, such as authors' names and related publications.

The [PDB File Format Documentaion](https://www.wwpdb.org/documentation/file-format) categorizes the data 
in a PDB file into twelve sections:

1. **Title**: Description of the experiment, biological molecules and their sources, and other general 
   information about the entry, such as identifiers, relationships to other entries, and bibliographic data.
2. **Remark**: Experimental details, annotations, comments, and information that do not match other predefined sections.
3. **Primary Structure**: Information on the primary structure of the polymers
   present in the PDB file, including the sequence of residues in each chain of the macromolecule(s);
   references to corresponding sequences in other databases, and possible conflicts between the two;
   and modifications of the residues in each sequence.
4. **Heterogen**: Description of non-standard groups (residues) in the entry, 
   such as prosthetic groups, inhibitors, solvent molecules, and ions, for which coordinates are supplied.
5. **Secondary Structure**: Information on the secondary structure elements, such as helices and sheets,
   found in the polymers present in the PDB file.
6. **Connectivity Annotation**: Description of disulfide bonds and other chemical connectivity data that is
   not implied by the primary structure.
7. **Miscellaneous Features**: Features within the macromolecule, such as an active site or environments
   surrounding a non-standard residue.
8. **Crystallographic**: Description of the geometry of the crystallographic experiment.
9. **Coordinate Transformation**: Coordinate system transformations.
10. **Coordinate**: The main piece of information detailing the atomic coordinates data.
11. **Connectivity**: 

12. **Bookkeeping**:


Record Format
-------------
Every PDB file is presented in a number of lines. Each line in the PDB entry file consists of 80
columns. The last character in each PDB entry should be an end-of- line indicator.
Each line in the PDB file is self-identifying. The first six columns of every line contains a record
name, that is left-justified and separated by a blank. The record name must be an exact match to
one of the stated record names in this format guide.
The PDB file may also be viewed as a collection of record types. Each record type consists of one
or more lines.
Each record type is further divided into fields.
For records that are fully described in fixed column format, columns not assigned to fields must be
left blank.

Types of Records
----------------
It is possible to group records into categories based upon how often the record type appears in an
entry.
One time, single line: There are records that may only appear one time and without continuations
in a file. It is an error for a duplicate of any of these records to appear in an entry.

One time, multiple lines: There are records that conceptually exist only once in an entry, but the
information content may exceed the number of columns available. These records are therefore
continued on subsequent lines. The second and subsequent lines contain a continuation field,
which is a right-justified integer.
This number increments by one for each additional line of the record, and is followed by a blank
character.

Multiple times, one line: Most record types appear multiple times, often in groups where the
information is not logically concatenated but is presented in the form of a list. Many of these record
types have a custom serialization that may be used not only to order the records, but also to
connect to other record types.

Multiple times, multiple lines: There are records that conceptually exist multiple times in an entry,
but the information content may exceed the number of columns available. These records are
therefore continued on subsequent lines.
The second and subsequent lines contain a continuation field which is a right-justified integer.
This number increments by one for each additional line of the record, and is followed by a blank
character.


Grouping: There are three record types used to group other records: MODEL, TER, ENDMDL
The MODEL/ENDMDL records surround groups of ATOM, HETATM, ANISOU, and TER records.
TER records indicate the end of a chain.

Other: The remaining record types have a detailed inner structure: JRNL, REMARK


Order of Records
----------------
All records in a PDB coordinate entry must appear in a defined order. Mandatory record types are
present in all entries. When mandatory data are not provided, the record name must appear in the
entry with a NULL indicator. Optional items become mandatory when certain conditions exist. Old
records that are not described here are deprecated.
* All records in a PDB coordinate entry must appear in a defined order. 
* Mandatory record types are present in all entries. 
* When mandatory data are not provided, the record name must appear in the entry with a NULL indicator. 
* Optional items become mandatory when certain conditions exist.



Basic Notions of the Format Description
---------------------------------------

Only non-control ASCII characters, as well as the space and end-of-line indicator, appear in a
PDB coordinate entry file. Namely:
abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ
1234567890
` - = [ ] \ ; ' , . / ~ ! @ # $ % ^ & * ( ) _ + { } | : " < > ?

Greek letters are spelled out, i.e., alpha, beta, gamma, etc.
Bullets are represented as (DOT).
Right arrow is represented as -->.
Left arrow is represented as <--.
If "=" is surrounded by at least one space on each side, then it is assumed to be an equal sign,
e.g., 2 + 4 = 6.
Commas, colons, and semi-colons are used as list delimiters in records that have one of the
following data types:
List
SList
Specification List
Specification

If a comma, colon, or semi-colon is used in any context other than as a delimiting character, then
the character must be escaped, i.e., immediately preceded by a backslash, "\".

Example - Use of “\” character:
COMPND MOL_ID: 1;
COMPND 2 MOLECULE: GLUTATHIONE SYNTHETASE;
COMPND 3 CHAIN: A;
COMPND 4 SYNONYM: GAMMA-L-GLUTAMYL-L-CYSTEINE\:GLYCINE LIGASE
COMPND 5 (ADP-FORMING);
COMPND 6 EC: 6.3.2.3;
COMPND 7 ENGINEERED: YES

COMPND MOL_ID: 1;
COMPND 2 MOLECULE: S-ADENOSYLMETHIONINE SYNTHETASE;
COMPND 3 CHAIN: A, B;
COMPND 4 SYNONYM: MAT, ATP\:L-METHIONINE S-ADENOSYLTRANSFERASE;
COMPND 5 EC: 2.5.1.6;
COMPND 6 ENGINEERED: YES;
COMPND 7 BIOLOGICAL_UNIT: TETRAMER;
COMPND 8 OTHER_DETAILS: TETRAGONAL MODIFICATION



### Contents
```{eval-rst}

+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Section                   | Record               | Description                                                                                                                                  | Mandatory                                                                                                           | Number of times | Number of lines |
+===========================+======================+==============================================================================================================================================+=====================================================================================================================+=================+=================+
| Title                     | HEADER               | First record of the entry, containing PDB ID, classification, and deposition date of the entry.                                              | Yes                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | OBSLTE               | Details about the removal of the entry from the database, and PDB IDs of entries that replaced it.                                           | When the entry has been replaced by a newer entry.                                                                  | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | TITLE                | Description of the experiment represented in the entry.                                                                                      | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SPLIT                | List of PDB IDs that compose a larger macromolecular complex, to which this entry belongs.                                                   | When the entry is split into multiple PDB entries.                                                                  | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | CAVEAT               | Severe error indicator.                                                                                                                      | When outstanding errors exists, e.g. regarding chirality.                                                           | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | COMPND               | Description of macromolecular contents of the entry.                                                                                         | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SOURCE               | Details about the biological source of macromolecules in the entry.                                                                          | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | KEYWDS               | List of keywords describing the macromolecule in the entry.                                                                                  | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | EXPDTA               | Experimental technique(s) used for structure determination.                                                                                  | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | NUMMDL               | Number of models in the entry.                                                                                                               | When more than one model exists, e.g. in NMR entries.                                                               | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | MDLTYP               | Additional annotations regarding the coordinates in the entry.                                                                               | When the entire polymer chain only contains backbone atoms, and in entries with NMR minimized average structures.   | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | AUTHOR               | List of contributors.                                                                                                                        | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REVDAT               | Release/revision history of the entry.                                                                                                       | Yes                                                                                                                 | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SPRSDE               | List of PDB IDs that were replaced (superseded) by this entry.                                                                               | When the entry has replaced other entries.                                                                          | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | JRNL                 | Literature citation that defines the entry.                                                                                                  | When a corresponding publication exists                                                                             | One             | Multiple        |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Rermark                   | REMARK 0             | Details of re-refinements performed on the structure, using data from an existing entry.                                                     | When the entry is a re-refined structure.                                                                           | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 1             | List of important publications related to the structure, chosen by the authors.                                                              | No                                                                                                                  | Multiple        | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 2             | Resolution of the model in diffraction studies.                                                                                              | Yes                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 3             | Information on final refinement programs used, and related statistics.                                                                       | Yes                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 4             | Data on the version of the PDB file format used to generate the file.                                                                        | N/A                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 5             | Obsolete statement, i.e. the reason for removal of the entry from the database.                                                              | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 100           | Data on deposition of the entry, including processing website and date.                                                                      | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 200           | Experimental details of X-ray diffraction studies.                                                                                           | When the structure was determined by single-crystal, fiber, or polycrystalline X-ray diffraction.                   | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 205           | Experimental details of fiber diffraction studies.                                                                                           | When the structure was determined by fiber diffraction, or non-crystalline sample studies.                          | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 210           | Experimental details of NMR studies.                                                                                                         | When the structure was determined by NMR studies.                                                                   | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 215           | Experimental details of NMR studies.                                                                                                         | When the structure was determined by solution NMR studies.                                                          | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 217           | Experimental details of NMR studies.                                                                                                         | When the structure was determined by solid-state NMR studies.                                                       | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 230           | Experimental details of neutron diffraction studies.                                                                                         | When the structure was determined by neutron diffraction studies.                                                   | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 240           | Experimental details of electron crystallography studies.                                                                                    | When the structure was determined by electron crystallography studies.                                              | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 245/247       | Experimental details of electron microscopy studies.                                                                                         | When the structure was determined by electron crystallography studies.                                              | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 250           | Experimental details of other types of studies.                                                                                              | When the structure was determined by studies other than those listed above (REMARKs 200–247).                       | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 265           | Experimental details of solution scattering studies.                                                                                         | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 280           | Information on the crystal, solvent content, Matthews coefficient, and crystallization conditions.                                           | When the structure was determined by single-crystal studies.                                                        | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 285           | Information on the unit cell of the crystal structure.                                                                                       | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 290           | Data on crystallographic symmetries.                                                                                                         | When the structure was determined by crystalline studies.                                                           | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 300           | Description of the biologically functional molecule.                                                                                         | When REMARK 350 exists.                                                                                             | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 350           | Description of all transformations needed to generate the biomolecule from the coordinates in the entry.                                     | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 375           | List of atoms which lie within 0.15 Å of a symmetry-related atom, and are thus considered to be on a special position.                       | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 400           | Further details about the macromolecular contents of the entry.                                                                              | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 450           | Further details about the biological source of the macromolecular contents of the entry.                                                     | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 465           | List of residues that are present in the SEQRES records but are completely absent from the coordinates section.                              | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 470           | List of non-hydrogen atoms of standard residues (or non-standard residues that are in SEQRES records) that are missing from the coordinates. | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 475           | List of residues that were modeled with zero occupancy.                                                                                      | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 480           | List of non-hydrogen atoms in residues modeled with zero occupancy.                                                                          | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 500           | Details about the stereochemistry of the structure.                                                                                          | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 525           | List of solvent atoms that are more than 5 Å away from any polymer chain.                                                                    | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 600           | Details on the heterogens in the entry.                                                                                                      | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 610           | List of non-polymer residues with missing atoms.                                                                                             | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 615           | List of non-polymer residues containing atoms with zero occupancy.                                                                           | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 620           | Details of metal coordination in the stucture.                                                                                               | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 630           | Description of (peptide) inhibitors that are present as a chemical component (HET group).                                                    | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 650           | Details on the helical portions of the entry.                                                                                                |                                                                                                                     | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 700           | Details on the sheet content of the structure.                                                                                               | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 800           | Details on important sites in the entry.                                                                                                     | When SITE records exist.                                                                                            | Multiple        | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 900           | Information on other related PDB entries.                                                                                                    | N/A                                                                                                                 | One             | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | REMARK 999           | Description of unusual/particular properties of polymer sequences in SEQRES records.                                                         | N/A                                                                                                                 | One             | Multiple        |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Primary Structure         | DBREF/DBREF1/DBREF2  | Reference to the entry in sequence databases.                                                                                                | When polymers exist (as opposed to only small molecules).                                                           | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SEQADV               | Identification of conflicts between the entry and the reference database.                                                                    | When sequence conflicts exist.                                                                                      | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SEQRES               | Primary sequence of backbone residues.                                                                                                       | When polymers exist.                                                                                                | Multiple        | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | MODRES               | Identification of modifications to standard residues.                                                                                        | When modified residues exist with coordinates data.                                                                 | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Heterogen                 | HET                  | Identification of non-standard (heterogen) groups.                                                                                           | When non-standard groups (other than water) exist.                                                                  | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | HETNAM               | Compound name of the heterogens.                                                                                                             | When non-standard groups (other than water) exist.                                                                  | Multiple        | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | HETSYN               | Synonymous compound names for heterogens.                                                                                                    | No                                                                                                                  | Multiple        | Multiple        |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | FORMUL               | Chemical formulas of non-standard groups.                                                                                                    | When non-standard groups (other than water) exist.                                                                  | Multiple        | Multiple        |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Secondary Structure       | HELIX                | Identification of helical substructures.                                                                                                     | No                                                                                                                  | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SHEET                | Identification of sheet substructures.                                                                                                       | No                                                                                                                  | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Connectivity Annotation   | SSBOND               | Identification of disulfide bonds.                                                                                                           | When disulfide bonds exist.                                                                                         | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | LINK                 | Identification of inter-residue bonds.                                                                                                       | When non-standard residues exist within a polymer.                                                                  | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | CISPEP               | Identification of peptide residues in cis conformation.                                                                                      | No                                                                                                                  | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Miscellaneous Features    | SITE                 | Identification of groups comprising important sites.                                                                                         | No                                                                                                                  | Multiple        | Multiple        |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Crystallographic          | CRYST1               | Crystallographic data, i.e. unit cell parameters, space group, and Z-value.                                                                  | Yes                                                                                                                 | One             | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Coordinate Transformation | ORIGX1/ORIGX2/ORIGX3 | Transformation from orthogonal coordinates to the submitted coordinates.                                                                     | Yes                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | SCALE1/SCALE2/SCLAE3 | Transformation from orthogonal coordinates to fractional crystallographic coordinates.                                                       | Yes                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | MTRIX1/MTRIX2/MTRIX3 | Transformations expressing non-crystallographic symmetry.                                                                                    | When the complete asymmetric unit must be generated from the given coordinates using non-crystallographic symmetry. | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Coordinate                | MODEL                | Specification of model number for multiple structures in a single coordinate entry.                                                          | When the entry contains more than one model.                                                                        | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | ATOM                 | Atomic coordinates for standard residues.                                                                                                    | When the entry contains polymers with standard residues.                                                            | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | ANISOU               | Anisotropic temperature factors.                                                                                                             | No                                                                                                                  | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | TER                  | Chain terminator.                                                                                                                            | When the entry contains polymers.                                                                                   | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | HETATM               | Atomic coordinates for non-standard (heterogen) groups.                                                                                      | When non-standard residues exist.                                                                                   | Multiple        | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | ENDMDL               | End-of-model indicator for multiple structures in a single coordinate entry.                                                                 | When MODEL records exist.                                                                                           | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Connectivity              | CONECT               | Connectivity records.                                                                                                                        | When either non-standard residues or disulfide bonds exist.                                                         | Multiple        | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
| Bookkeeping               | MASTER               | Control record for bookkeeping and validation; contains checksums (i.e. counts) of other records.                                            | Yes                                                                                                                 | One             | One             |
|                           +----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+
|                           | END                  | End-of-file identifier.                                                                                                                      | Yes                                                                                                                 | One             | One             |
+---------------------------+----------------------+----------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------+-----------------+-----------------+

```

```{bibliography}
:style: unsrt
```