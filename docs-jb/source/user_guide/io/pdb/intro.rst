Overview
===================================


The Protein Data Bank (PDB) is an archive of experimentally determined three-dimensional
structures of biological macromolecules that serves a global community of researchers, educators,
and students. The data contained in the archive include atomic coordinates, crystallographic
structure factors and NMR experimental data. Aside from coordinates, each deposition also
includes the names of molecules, primary and secondary structure information, sequence
database references, where appropriate, and ligand and biological assembly information, details
about data collection and structure solution, and bibliographic citations.

This comprehensive guide describes the "PDB format" used by the members of the worldwide
Protein Data Bank (wwPDB; Berman, H.M., Henrick, K. and Nakamura, H. Announcing the
worldwide Protein Data Bank. Nat Struct Biol 10, 980 (2003)).

Version 3.30: Current version, minor addenda to Version 3.2. , introducing a small number of
changes and extensions supporting the annotation practices adopted by the wwPDB beginning in
July 13 2011 including REMARK 0, REMARK 3, REMARK 400 and REMARK 630.
July 13 2011, initial version 3.30.
October 05, 2011, update REMARK 350 and OBSOLETE.
May 09, 2012, update description in REMARK 470, HET, ATOM and MODEL.
November 21, 2012, minor typo corrections.


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

.. include:: pdb_data_table.rst