"""

"""

from typing import Optional

import numpy as np
import pandas as pd

from opencadd.io.pdb import records, fields


class Title:
    """
    This section contains records used to describe the experiment and the biological macromolecules
    present in the entry: HEADER, OBSLTE, TITLE, SPLIT, CAVEAT, COMPND, SOURCE,
    KEYWDS, EXPDTA, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and JRNL records.

    References
    ----------
    https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html
    """
    def __init__(
            self,
            header: records.Header,
            title: str,
            compound: pd.DataFrame,

    ):
        self.header: records.Header = header
        self.title: str = title
        self.compound: pd.DataFrame = compound
        return

    @property
    def header(self) -> records.Header:
        return self._header

    @header.setter
    def header(self, value):
        if not isinstance(value, records.Header):
            raise TypeError(
                "Attribute `header` expects a `Header`. "
                f"Input value {value} had type {type(value)}."
            )
        self._header = value
        return

    @property
    def title(self) -> str:
        """
        Title of the experiment.
        It is a free text that describes the contents of the entry and any procedures or
        conditions that distinguish it from similar entries. Some data that may be included are
        e.g. experiment type, description of the mutation, and the fact that only alpha carbon
        coordinates have been provided in the entry.
        """
        return self._title

    @title.setter
    def title(self, value):
        if not isinstance(value, str):
            raise TypeError(
                "Attribute `title` expects a string. "
                f"Input value {value} had type {type(value)}."
            )
        self._title = value
        return

    @property
    def compound(self) -> pd.DataFrame:
        """
        The COMPND record describes the macromolecular contents of an entry. Some cases where the
        entry contains a standalone drug or inhibitor, the name of the non-polymeric molecule will
        appear in this record.
        """
        return self._compound

    @compound.setter
    def compound(self, value):
        if not isinstance(value, pd.DataFrame):
            raise TypeError(
                "Attribute `compound` expects a `pandas.DataFrame`. "
                f"Input value {value} had type {type(value)}."
            )
        self._compound = value
        return

    @property
    def pdb_format(self) -> str:
        return

    @property
    def pdb_format_title(self) -> str:
        """
        TITLE Record formatted as in a PDB file.
        """
        num_lines_needed = int(np.ceil(len(self.title) / 70))
        continuation = fields.Continuation.to_pdb(num_lines=num_lines_needed)
        title = [self.title[i:i+70] for i in range(0, len(self._title), 70)]
        lines = [
            f"TITLE{'':3}{continuation}{title:<70}"
            for continuation, title in zip(continuation, title)
        ]
        return "\n".join(lines)

    @property
    def pdb_format_compound(self) -> str:
        spec_str = ""
        for mol_id in self.compound.index:
            mol_data = self.compound.loc[mol_id].values
            mask_nan = np.logical_not(np.isnan(mol_data))
            non_nan_data = mol_data[mask_nan]
            spec_str += f"MOL_ID: {mol_id};"
            for idx, token in enumerate(self.compound.column.values[mask_nan]):
                if token in ("CHAIN", "SYNONYM", "EC"):
                    spec_str += f"{token}: {', '.join(non_nan_data[idx])};"
                elif token in ("ENGINEERED", "MUTATION"):
                    spec_str += f"{token}: {'YES' if non_nan_data[idx] else 'NO'};"
                else:
                    spec_str += f"{token}: {non_nan_data[idx]};"

        lines_needed = int(np.ceil(len(spec_str) / 70))
        continuations = ["  "] + [f"{line_num:>3}" for line_num in range(2, lines_needed + 1)]
        specs = [spec_str[i:i + 70] for i in range(0, len(spec_str), 70)]
        return "\n".join(
            [
                f"COMPND{'':1}{continuation}{title:<70}"
                for continuation, title in zip(continuations, specs)
            ]
        )


class Remark:
    pass


class PrimaryStructure:
    """
    Information on the primary structure of the polymers present in the PDB file.
    Contains the sequence of residues in each chain of the macromolecule(s). Embedded in these
    records are chain identifiers and sequence numbers that allow other records to link into the sequence.
    Peptide and/or nucleotide sequence and the relationship between the PDB sequence and that
    found in the sequence database(s).
    """
    def __init__(
            self,
            dbref: Optional[pd.DataFrame] = None,
            seqadv: Optional[pd.DataFrame] = None,
            seqres: Optional[pd.DataFrame] = None,
            modres: Optional[pd.DataFrame] = None,
    ):
        self._dbref: pd.DataFrame = dbref
        self._seqadv: pd.DataFrame = seqadv
        self._seqres: pd.DataFrame = seqres
        self._modres: pd.DataFrame = modres
        return

    @property
    def dbref(self) -> Optional[pd.DataFrame]:
        """
        The DBREF and DBREF1/DBREF2 records of the PDB file.
        Provides cross-references between the sequences of the polymers (chains)
        in the PDB file (as it appears in the SEQRES record), and a corresponding database
        sequence.

        Returns
        -------
        pandas.DataFrame | None
            The columns are defined as follows:
            chain_id: Chain identifier of the polymer in the PDB file.
            residue_num_begin: Initial residue sequence number of the polymer in the PDB file.
            residue_icode_begin: Initial residue insertion code of the polymer in the PDB file.
            residue_num_end: Ending residue sequence number of the polymer in the PDB file.
            residue_icode_end: Ending residue insertion code of the polymer in the PDB file.
            db: Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            db_chain_accession: Accession code of the polymer in the database.
            db_chain_id: Reference to 'chain_id' in the database.
            db_residue_num_begin: Reference to 'residue_num_begin' in the database.
            db_residue_icode_begin: Reference to 'residue_icode_begin' in the database.
            db_residue_num_end: Reference to 'residue_num_end' in the database.
            db_residue_icode_end: Reference to 'residue_icode_end' in the database.

        Notes
        -----
        Since both DBREF and DBREF1/DBREF2 records contain the same type of information,
        they are both parsed into the same dataframe.
        """
        return self._dbref

    @property
    def seqadv(self) -> Optional[pd.DataFrame]:
        """
        The SEQADV record of the PDB file.
        Identifies the differences between sequence information in the SEQRES records of the PDB
        entry and the sequence database entry given in DBREF.
        No assumption is made as to which database contains the correct data.

        Returns
        -------
        pandas.DataFrame | None
            The columns are defined as follows:
            residue_name: Name of the conflicting residue in the PDB file.
            chain_id: Chain identifier of the conflicting residue's parent polymer in the PDB file.
            residue_num: Sequence number of the conflicting residue in the PDB file.
            residue_icode: Insertion code of the conflicting residue in the PDB file.
            db: Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            db_chain_accession: Accession code of the polymer (chain) in the database.
            db_residue_name: Reference to 'residue_name' in the database.
            db_residue_num: Reference to 'residue_num' in the database.
            description: Description of the conflict. Some possible comments are:
                'Cloning artifact', 'Expression tag', 'Conflict', 'Engineered', 'Variant',
                'Insertion', 'Deletion', 'Microheterogeneity', 'Chromophore'.
                If a conflict is not classifiable by these terms, a reference to either a
                published paper, a PDB entry, or a REMARK within the entry is given.
                The comment 'SEE REMARK 999' is used when the comment is too long.
        """
        return self._seqadv

    @property
    def seqres(self) -> Optional[pd.DataFrame]:
        """
        The SEQRES records of the PDB file.
        Contains information on the sequence of each polymer (i.e. chain) in the file.
        The components (i.e. residues) of each sequence may be standard or modified amino/nucleic
        acids, or other residues that are linked to the standard backbone in the polymer.
        Components that are linked to side-chains, or sugars and/or bases are not listed here.
        Residues in the ATOM records must agree with the corresponding sequence in SEQRES records.

        Returns
        -------
        pandas.DataFrame | None
            The columns are defined as follows:
            chain_id: Chain identifier of the polymer in the PDB file.
            residue_count: Number of residues in the polymer (chain).
            residue_names: Name of the residues in the polymer (chain).

        Notes
        -----
        * Ribo- and deoxyribonucleotides are distinguished; ribo residues are identified with the
        residue names A, C, G, U and I, while deoxy residues are identified with the residue names
        DA, DC, DG, DT and DI. Modified nucleotides are marked by separate 3-letter residue codes.
        * Known problems:
          * Polysaccharides are not properly represented.
          * If the starting position of a sequence is unknown, the sequence cannot be described.
          * For cyclic peptides, a random residue must be assigned as the N-terminus.
        """
        return self._seqres

    @property
    def modres(self) -> Optional[pd.DataFrame]:
        """
        The MODRES records of the PDB file.
        Provides descriptions of modifications (e.g. chemical or post-translational) to protein
        and nucleic acid residues. Included are correlations between residue names given in a
        PDB entry and standard residues.
        Residues modified post-translationally, enzymatically, or by design are described.
        In those cases where the wwPDB has opted to use a non-standard residue name for the
        residue, MODRES also correlates the new name to the precursor standard residue name.

        Returns
        -------
        pandas.DataFrame | None
            The columns are defined as follows:
            residue_name: Name of the modified residue, as used in the PDB file.
            chain_id: Chain identifier of the modified residue's parent chain in the PDB file.
            residue_num: Sequence number of the modified residue in the PDB file.
            residue_icode: Insertion code of the modified residue in the PDB file.
            residue_name_std: Standard name of the modified residue.
            description: Description of the modification.

        Notes
        -----
        * D-amino acids are given their own residue name (resName), i.e., DAL for D-alanine. The
        residue name appears in the SEQRES records, and has the associated MODRES, HET, and FORMUL
        records. The coordinates are given as HETATMs within the ATOM records and occur in the
        correct order within the chain. This ordering is an exception to the Order of Records.

        * When a standard residue name is used to describe a modified site, residue_name
        residue_name_std contain the same value.

        * MODRES is mandatory when modified standard residues exist in the entry, but is not
        required if coordinate records are not provided for the modified residue.
        """
        return self._modres


class Heterogen:
    """
    The heterogen section of a PDB formatted file contains the complete description of non-standard
    residues in the entry. Detailed chemical definitions of non-polymer chemical components are
    described in the Chemical Component Dictionary (ftp://ftp.wwpdb.org/pub/pdb/data/monomers)
    """
    def __init__(
            self,
            heterogens: pd.DataFrame,
            occurrences: pd.DataFrame
    ):
        self._heterogen = heterogens
        self._occurrences = occurrences
        return

    @property
    def het(self) -> pd.DataFrame:
        """
        HET records are used to describe non-standard residues, such as prosthetic groups, inhibitors,
        solvent molecules, and ions for which coordinates are supplied. Groups are considered HET if
        they are not part of a biological polymer described in SEQRES and considered to be a molecule
        bound to the polymer, or they are a chemical species that constitute part of a biological polymer
        and is not one of the following:
        • standard amino acids, or
        • standard nucleic acids (C, G, A, U, I, DC, DG, DA, DU, DT and DI), or
        • unknown amino acid (UNK) or nucleic acid (N) where UNK and N are used to indicate the
        unknown residue name.

        HET records also describe chemical components for which the chemical identity is unknown, in
        which case the group is assigned the hetID UNL (Unknown Ligand).

        Returns
        -------
        pandas.DataFrame
            heterogen_id: Identifier of the non-standard residue; each unique ID represents a
                unique molecule.
            chain_id: Chain identifier of the non-standard residue's parent chain in the PDB file.
            residue_num: Sequence number of the non-standard residue in the PDB file.
            residue_icode: Insertion code of the non-standard residue in the PDB file.
            hetatm_count: Number of HETATM records present in the PDB file corresponding to this
                molecule.
            description: Description of the non-standard residue.

        Notes
        -----
        * There is a separate HET record for each occurrence of the HET group in an entry.
        * A particular HET group is represented in the PDB archive with a unique hetID.
        * PDB entries do not have HET records for water molecules, deuterated water, or methanol (when
        used as solvent).
        * Unknown atoms or ions will be represented as UNX with the chemical formula X1. Unknown
        ligands are UNL; unknown amino acids are UNK.
        """


class SecondaryStructure:
    pass


class Connectivity:
    pass


class Crystallographic:
    pass



class Coordinate:

    def __init__(
            self,
            atoms_data: pd.DataFrame,
            ter
    ):
        pass


class Bookkeeping:
    """
    Bookkeeping section of a PDB file. It contains the information from the MASTER record,
    which is a control record, counting the occurrence (checksum) of selected record types
    in the PDB file.

    Attributes
    ----------
    remark : int
        Number of REMARK records.
    het : int
        Number of HET records.
    helix : int
        Number of HELIX records.
    sheet : int
        Number of SHEET records.
    site : int
        Number of SITE records.
    xform : int
        Number of coordinate transformation records, i.e. ORIGX, SCALE and MTRIX records combined.
    coord : int
        Number of atomic coordinate records, i.e. ATOM and HETATM records combined.
    ter : int
        Number of TER records.
    connect : int
        Number of CONECT records.
    seqres : int
        Number of SEQRES records.

    Notes
    -----
    When there are multiple models in the file, only the first model is counted.
    """
    remark: int
    het: int
    helix: int
    sheet: int
    site: int
    xform: int
    coord: int
    ter: int
    connect: int
    seqres: int