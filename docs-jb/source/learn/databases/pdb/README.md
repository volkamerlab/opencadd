## PDB File Description

### Reference
* [Protein Data Bank Contents Guide: Atomic Coordinate Entry Format Description, Version 3.30.](https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf)

### Characters

* Allowed are only non-control ASCII characters, plus space and end-of-line indicator, i.e.:
  * abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ 1234567890`-=[]\;',./~!@#$%^&*()_+{}|:"<>?
* Commas, colons, and semicolons are only used as list delimeters. In any other context, they 
  must be escaped, i.e. immediately preceded by a backslash, e.g. \\;

### Structure
* Every line is 80 characters (followed by an end-of-line indicator).
* First 6 characters of every line contain the *record* name (left-justified).


* Records can be grouped based on their frequency and number of lines:
  #### One time, single line: 
    * Only appear once in the document, on a single line.
    * It is thus an error for a duplicate to appear.
    * These are:
      * **CRYST1**: Unit cell parameters, space group, and Z.
      * **END**: Last record in the file. 
      * **HEADER**: First line of the entry, contains PDB ID code, classification, and date of 
        deposition. 
      * **NUMMDL**: Number of models. 
      * **MASTER**: Control record for bookkeeping. 
      * **ORIGXn** (n = 1, 2, or 3): Transformation from orthogonal coordinates to the submitted 
        coordinates. 
      * **SCALEn** (n = 1, 2, or 3): Transformation from orthogonal coordinates to fractional 
        crystallographic 
        coordinates.

  #### One time, multiple lines:
    * Only appear once in the document, but can span multiple lines.
    * All lines other than the first line start with a right-justified integer, which 
      starts at 2 and increments by one for each additional line of the record, and is followed 
      by a blank character.
    * These are:
      * **AUTHOR**: List of contributors. 
      * **CAVEAT**: Severe error indicator. 
      * **COMPND**: Description of macromolecular contents of the entry. 
      * **EXPDTA**: Experimental technique used for the structure determination. 
      * **MDLTYP**: Contains additional annotation pertinent to the coordinates presented in the 
        entry. 
      * **KEYWDS**: List of keywords describing the macromolecule. 
      * **OBSLTE**: Statement that the entry has been removed from distribution and list of the ID 
        code(s) which replaced it. 
      * **SOURCE**: Biological source of macromolecules in the entry. 
      * **SPLIT**: List of PDB entries that compose a larger macromolecular complexes. 
      * **SPRSDE**: List of entries obsoleted from public release and replaced by current entry.
      * **TITLE**: Description of the experiment represented in the entry.

  #### Multiple times, one line
    * Appear multiple times (often in groups), each time only on a single line.
    * Represent tabular data, and comprise the most of the file.
    * These are:
      * **ANISOU**: Anisotropic temperature factors. 
      * **ATOM**: Atomic coordinate records for standard groups. 
      * **CISPEP**: Identification of peptide residues in cis conformation. 
      * **CONECT**: Connectivity records. 
      * **DBREF**: Reference to the entry in the sequence database(s). 
      * **HELIX**: Identification of helical substructures. 
      * **HET**: Identification of non-standard groups heterogens). 
      * **HETATM**: Atomic coordinate records for heterogens. 
      * **LINK**: Identification of inter-residue bonds. 
      * **MODRES**: Identification of modifications to standard residues. 
      * **MTRIXn** (n = 1, 2, or 3): Transformations expressing non-crystallographic symmetry. 
        There may be multiple sets of these records. 
      * **REVDAT**: Revision date and related information. 
      * **SEQADV**: Identification of conflicts between PDB and the named sequence database. 
      * **SHEET**: Identification of sheet substructures. 
      * **SSBOND**: Identification of disulfide bonds.

  #### Multiple times, multiple lines:
    * Appear multiple times, and each time the information may span multiple lines.
    * Just like in 'one time, multiple lines', all lines other than the first line start with a 
      right-justified integer, which starts at 2 and increments by one for each additional line of the record, and is followed 
      by a blank character.
  
  #### Grouping
    * There are three record types used to group other records:
      * **MODEL**: Specification of model number for multiple structures in a single coordinate entry.
      * **ENDMDL**: End-of-model record for multiple structures in a single coordinate entry.
      * **TER**: Chain terminator (indicates the end of a chain).
    * MODEL/ENDMDL records surround groups of ATOM, HETATM, ANISOU, and TER records.
  
  #### Structured
    * There are two record types with a detailed inner structure:
      * **JRNL**: Literature citation that defines the coordinate set. 
      * **REMARK**: General remarks; they can be structured or free form.

 

