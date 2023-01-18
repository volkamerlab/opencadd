"""

"""
import datetime
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from opencadd.io.pdb import fields


class Title:
    """
    Title section of the PDB file, describing the experiment, biological molecules, and other general
    information about the entry, corresponding publications and authors.

    Notes
    -----
    * This section contains all data from the PDB records HEADER, OBSLTE, TITLE, SPLIT, CAVEAT, COMPND, SOURCE,
    KEYWDS, EXPDTA, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and JRNL.

    References
    ----------
    https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html
    """
    def __init__(
            self,
            pdb_id: Optional[str] = None,
            deposition_date: Optional[datetime.date] = None,
            classification: Optional[Tuple[Tuple[str]]] = None,
            replacement_date: Optional[datetime.date] = None,
            replacement_pdb_ids: Optional[Sequence[str]] = None,
            title: Optional[str] = None,
            split_pdb_ids: Optional[Sequence[str]] = None,
            caveat: Optional[str] = None,
            compound: Optional[pd.DataFrame] = None,
            source: Optional[pd.DataFrame] = None,
            keywords: Optional[Sequence[str]] = None,
            experimental_techniques: Optional[Sequence[str]] = None,
            count_models: Optional[int] = None,
            structural_annotations: Optional[Sequence[str]] = None,
            authors: Optional[Sequence[str]] = None,
    ):
        self._pdb_id: str = pdb_id
        self._deposition_date: datetime.date = deposition_date
        self._classification: Tuple[Tuple[str]] = classification
        self._replacement_date: datetime.date = replacement_date
        self._replacement_pdb_ids: Tuple[str] = tuple(replacement_pdb_ids)
        self._title: str = title
        self._split_pdb_ids: Tuple[str] = tuple(split_pdb_ids)
        self._caveat: str = caveat
        self._compound: pd.DataFrame = compound
        self._source: pd.DataFrame = source
        self._keywords: Tuple[str] = tuple(keywords)
        self._experimental_techniques: Tuple[str] = tuple(experimental_techniques)
        self._count_models: int = count_models
        self._structural_annotations: Tuple[str] = tuple(structural_annotations)
        self._authors: tuple[str] = tuple(authors)
        return

    def __repr__(self):
        return f"Header({self.pdb_id}, {self.deposition_date}, {self.classification})"

    def __str__(self):
        lines = [
            f"PDB-ID: {self.pdb_id}",
            f"Deposition Date: {self.deposition_date}",
            "Classification (per entity):"
        ] + [f"\t{i+1}. {', '.join(entity)}" for i, entity in enumerate(self.classification)]
        return "\n".join(lines)

    @property
    def pdb_id(self) -> str:
        """
        PDB identification code (PDB ID) of the entry.
        This identifier is unique within the Protein Data Bank (PDB).

        Returns
        -------
        str | None
            The PDB identification code (PDB ID), consisting of 4 characters, the first of which is a digit
            in the range 0 - 9; the remaining 3 are alphanumeric, and letters are upper case only.

        Notes
        -----
        * PDB IDs with a 0 as the first character do not contain coordinate data; these so-called
        "no-coordinates" (NOC) files were removed from the PDB archive.
        * This property corresponds to the 'idCode' field of the HEADER record in the PDB file.
        """
        return self._pdb_id

    @property
    def deposition_date(self) -> datetime.date:
        """
        Date of deposition of the entry at the Protein Data Bank.

        Returns
        -------
        datetime.date | None

        Notes
        -----
        * This property corresponds to the 'depDate' field of the HEADER record in the PDB file.
        """
        return self._deposition_date

    @property
    def classification(self) -> Tuple[Tuple[str]]:
        """
        Classification of each molecule within the entry.

        Returns
        -------
        Tuple[Tuple[str]] | None
            Each sub-tuple corresponds to one molecule in the entry, with each element
            describing one classification/function.

        Notes
        -----
        * Classification may be based on function, metabolic role, molecule type, cellular location
        etc., but it exactly matches one of a collection of strings, available at:
        https://www.wwpdb.org/docs.html
        *  Due to the limited length of the classification field, strings must sometimes be
        abbreviated. In these cases, the full terms are given in KEYWDS records in no strict order.
        * This property corresponds to the 'classification' field of the HEADER record in the PDB file.
        """
        return self._classification

    @property
    def replacement_date(self) -> Optional[datetime.date]:
        """
        The date the entry was removed (“obsoleted”) from the PDB's full release.
        This only appears in entries that have been removed from public distribution.

        Returns
        -------
        datetime.date | None

        Notes
        -----
        * All OBSLTE entries are available from the PDB archive at:
        https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete.
        * This property corresponds to the 'repDate' field of the OBSLTE records in the PDB file.
        """
        return self._replacement_date

    @property
    def replacement_pdb_ids(self) -> Optional[np.ndarray]:
        """
        PDB IDs of the new entries that have replaced this entry.
        This only appears in entries that have been removed from public distribution.

        Returns
        -------
        numpy.ndarray[ndims:1, dtype:<U4] | None

        Notes
        -----
        * All OBSLTE entries are available from the PDB archive at:
        https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete.
        * This property corresponds to the 'rIdCode' fields of the OBSLTE records in the PDB file.
        """
        return self._replacement_pdb_ids

    @property
    def title(self) -> str:
        """
        Title of the experiment, describing the contents of the entry and any procedures or
        conditions that distinguish it from similar entries.
        Some data that may be included are e.g. experiment type, description of the mutation,
        and the fact that only alpha carbon coordinates have been provided in the entry.

        Returns
        -------
        str | None

        Notes
        -----
        * This property corresponds to the 'title' field of the TITLE record in the PDB file.
        """
        return self._title

    @property
    def split_pdb_ids(self) -> Optional[np.ndarray]:
        """
        PDB IDs that are required to reconstitute a complete complex.
        This only appears in instances where the entry composes a part of a large macromolecular complex.

        Returns
        -------
        numpy.ndarray[ndims:1, dtype:<U4] | None

        Notes
        -----
        * If this property exist, REMARK 350 will contain an amended statement to reflect the entire complex.
        * This property corresponds to the 'idCode' fields of the SPLIT records in the PDB file.
        """
        return self._split_pdb_ids

    @property
    def caveat(self) -> str:
        """
        A free text comment that warns of errors and unresolved issues in the entry, if any.
        It will also be included in cases where the PDB is unable to verify the transformation of the
        coordinates back to the crystallographic cell. In these cases, the molecular structure may still
        be correct.

        Returns
        -------
        str | None

        Notes
        -----
        * This property corresponds to the 'comment' field of the CAVEAT records in the PDB file.
        """
        return self._caveat

    @property
    def compound(self) -> pd.DataFrame:
        """
        Description of the macromolecular contents of the PDB file, or a standalone drug or inhibitor in cases
        where the entry does not contain a polymer.

        Returns
        -------
        pandas.DataFrame
            Index:
            * mol_id (int): Enumerates each molecule; the same ID appears also in the `source` property of
            this object (i.e. `Title.source`).
            Columns:
            * name (str): Name of the (macro)molecule. For chimeric proteins, the protein name is
            comma-separated and may refer to the presence of a linker, e.g. "protein_1, linker, protein_2".
            * chain_ids (numpy.ndarray, dtype: <U1): Chain identifiers in the macromolecule.
            * fragment (str): Name or description of a domain or region of the molecule.
            * synonyms (numpy.ndarray, dtype: str): synonyms for the molecule's name.
            * enzyme_commission_num (numpy.ndarray, dtype: str): Enzyme commision (EC) numbers associated with
            the molecule.
            * engineered (bool): Whether the molecule was produced using recombinant technology or by purely
            chemical synthesis.
            * mutation (bool): Whether there is a mutation in the molecule.
            * description (str): Additional free-text comment.

        Notes
        -----
        * For one (macro)molecule, multiple entries may exist in the dataframe, where each entry corresponds
        to a certain 'fragment' inside the molecule.
        * For nucleic acids, 'name' may contain asterisks, which are for ease of reading.
        * When residues with insertion codes occur in 'fragment' and 'description' the insertion code must be
        given in square brackets, e.g. "H57[A]N".
        * This property corresponds to the 'compound' fields of the COMPND records in the PDB file.
        The 'compound' field is a specification list, with a defined set of tokens for each component. These
        tokens correspond to the columns (or the index) of the returned dataframe, as follows:
            * MOL_ID: mol_id (index)
            * MOLECULE: name
            * CHAIN: chain_ids
            * FRAGMENT: fragment
            * SYNONYM: synonyms
            * EC: enzyme_commission_num
            * ENGINEERED: engineered
            * MUTATION: mutation
            * OTHER_DETAILS: description
        """
        return self._compound

    @property
    def source(self) -> pd.DataFrame:
        """
        Information on the biological/chemical source of each biological molecule in the PDB file, or
        a standalone drug or inhibitor in cases where the entry does not contain a polymer.

        Returns
        -------
        pandas.DataFrame | None
            Index:
            * mol_id: Enumerates each molecule; the same ID appears also in the `compound` property of
            this object (i.e. `Title.compound`).
            Columns:
            * synthetic: Indicates a chemically synthesized source.
            * fragment: Specifies a domain or fragment of the molecule.
            * organism: Common name of the organism
            * organism_sci: Scientific name of the organism.
            * organism_tax_id: NCBI Taxonomy ID of the organism.
            * strain: Identifies the strain.
            * variant: Identifies the variant.
            * cell_line: The specific line of cells used in the experiment.
            * atcc_id: American Type Culture Collection tissue culture number.
            * organ: Organized group of tissues that carries on a specialized function.
            * tissue: Organized group of cells with a common function and structure.
            * cell: Identifies the particular cell type.
            * organelle: Organized structure within a cell.
            * secretion: Identifies the secretion, such as saliva, urine, or venom, from which the molecule
            was isolated.
            * cell_loc: Identifies the location inside/outside the cell, where the compound was found.
            Examples are: 'extracellular', 'periplasmic', 'cytosol'.
            * plasmid: Identifies the plasmid containing the gene.
            * gene: Identifies the gene.
            * expsys: Expression system, i.e. common name of the organism in which the molecule was expressed.
            * expsys_sci: Scientific name of the expression system.
            * expsys_tax_id: NCBI Taxonomy ID of the expression system.
            * expsys_strain: Strain of the organism in which the molecule was expressed.
            * expsys_variant: Variant of the organism used as the expression system.
            * expsys_cell_line: The specific line of cells used as the expression system.
            * expsys_atcc_id: American Type Culture Collection tissue culture number of the expression system.
            * expsys_organ: Specific organ which expressed the molecule.
            * expsys_tissue: Specific tissue which expressed the molecule.
            * expsys_cell: Specific cell type which expressed the molecule.
            * expsys_organelle: Specific organelle which expressed the molecule.
            * expsys_cell_loc: Identifies the location inside or outside the cell which expressed the molecule.
            * expsys_vector_type: Identifies the type of vector used, i.e. plasmid, virus, or cosmid.
            * expsys_vector: Identifies the vector used.
            * expsys_plasmid: Plasmid used in the recombinant experiment.
            * expsys_gene: Name of the gene used in recombinant experiment.
            * details: Other details about the source.

        Notes
        -----
        * Sources are described by both the common name and the scientific name, e.g., genus and species.
        Strain and/or cell-line for immortalized cells are given when they help to uniquely identify
        the biological entity studied.
        * Molecules prepared by purely chemical synthetic methods are identified by the
        column `synthetic` with a "YES" value, or an optional value, such as "NON-BIOLOGICAL
        SOURCE" or "BASED ON THE NATURAL SEQUENCE". The `engineered` column in the COMPND record
        is also set in such cases.
        * Hybrid molecules prepared by fusion of genes are treated as multi-molecular systems for
        the purpose of specifying the source. The column `fragment` is used to associate the source
        with its corresponding fragment.
        * When necessary to fully describe hybrid molecules, tokens may appear more than once for
        a given `mol_id`.
        * This property corresponds to the 'srcName' fields of the SOURCE records in the PDB file.
        The 'srcName' field is a specification list, with a defined set of tokens for each component. These
        tokens correspond to the columns (or the index) of the returned dataframe, as follows:
            * MOL_ID: mol_id (index)
            * SYNTHETIC: synthetic
            * FRAGMENT: fragment
            * ORGANISM_COMMON: organism
            * ORGANISM_SCIENTIFIC: organism_sci
            * ORGANISM_TAXID: organism_tax_id
            * STRAIN: strain
            * VARIANT: variant
            * CELL_LINE: cell_line
            * ATCC: atcc_id
            * ORGAN: organ
            * TISSUE: tissue
            * CELL: cell
            * ORGANELLE: organelle
            * SECRETION: secretion
            * CELLULAR_LOCATION: cell_loc
            * PLASMID: plasmid
            * GENE: gene
            * EXPRESSION_SYSTEM_COMMON: expsys
            * EXPRESSION_SYSTEM: expsys_sci
            * EXPRESSION_SYSTEM_TAXID: expsys_tax_id
            * EXPRESSION_SYSTEM_STRAIN: expsys_strain
            * EXPRESSION_SYSTEM_VARIANT: expsys_variant
            * EXPRESSION_SYSTEM_CELL_LINE: expsys_cell_line
            * EXPRESSION_SYSTEM_ATCC_NUMBER: expsys_atcc_id
            * EXPRESSION_SYSTEM_ORGAN: expsys_organ
            * EXPRESSION_SYSTEM_TISSUE: expsys_tissue
            * EXPRESSION_SYSTEM_CELL: expsys_cell
            * EXPRESSION_SYSTEM_ORGANELLE: expsys_organelle
            * EXPRESSION_SYSTEM_CELLULAR_LOCATION: expsys_cell_loc
            * EXPRESSION_SYSTEM_VECTOR_TYPE: expsys_vector_type
            * EXPRESSION_SYSTEM_VECTOR: expsys_vector
            * EXPRESSION_SYSTEM_PLASMID: expsys_plasmid
            * EXPRESSION_SYSTEM_GENE: expsys_gene
            * OTHER_DETAILS: details
        """
        return self._source

    @property
    def keywords(self) -> np.ndarray:
        """
        Keywords/terms relevant to the PDB file, similar to that found in journal articles.
        The provided terms may for example describe functional classification, metabolic role, known
        biological or chemical activity, or structural classification.

        Returns
        -------
        numpy.ndarray[ndims: 1, dtype: str]
            Array of keywords.

        Notes
        -----
        * The classifications given in the `classification` property of this object
        (i.e. `Title.classification`) are also repeated here, with two differences: Unlike in `classification`,
        here the keywords are not grouped per molecule, but they are given unabbreviated.
        * This property corresponds to the 'keywds' fields of the KEYWDS records in the PDB file.
        """
        return self._keywords

    @property
    def experimental_techniques(self) -> np.ndarray:
        """
        Experimental techniques used to obtain the structure.

        Returns
        -------
        numpy.ndarray[ndims: 1, dtype: str]
            Array of experimental techniques. Permitted values include:
            "X-RAY DIFFRACTION", "FIBER DIFFRACTION", "NEUTRON DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY",
            "ELECTRON MICROSCOPY", "SOLID-STATE NMR", "SOLUTION NMR" and "SOLUTION SCATTERING".

        Notes
        -----
        * Since October 15, 2006, theoretical models are no longer accepted for deposition. Any
        theoretical models deposited prior to this date are archived at:
        ftp://ftp.wwpdb.org/pub/pdb/data/structures/models
        * This property corresponds to the 'technique' fields of the EXPDATA records in the PDB file.
        """
        return self._experimental_techniques

    @property
    def count_models(self) -> int:
        """
        Total number of models in the PDB file.

        Returns
        -------
        int | None

        Notes
        -----
        * This property corresponds to the 'modelNumber' field of the NUMMDL record in the PDB file.
        """
        return self._count_models

    @property
    def structural_annotations(self) -> Tuple[str]:
        """
        Additional structural annotations pertinent to the coordinates in the PDB file, used to highlight
        certain features.

        Returns
        -------
        tuple[str]

        Notes
        -----
        *  For entries that are determined by NMR methods and the coordinates deposited are either a
        minimized average or regularized mean structure, the tag "MINIMIZED AVERAGE" will be present as the
        first element of the returned tuple.
        * Where the entry contains entire polymer chains that have only either C-alpha (for proteins) or
        P atoms (for nucleotides), the contents of such chains will be described along with the
        chain identifier, e.g. " CA ATOMS ONLY, CHAIN A, B". For these polymeric chains,
        REMARK 470 (Missing Atoms) will be omitted.
        * This property corresponds to the 'comment' fields of the MDLTYP record in the PDB file.
        """
        return self._structural_annotations

    @property
    def authors(self) -> Tuple[str]:
        """
        Names of the people responsible for the contents of the PDB file.

        Returns
        -------
        tuple[str]

        Notes
        -----
        * First and middle names are indicated by initials, each followed by a period, and precede the surname.
        * Only the surname (family or last name) of the author is given in full.
        * Hyphens can be used if they are part of the author's name.
        * Apostrophes are allowed in surnames.
        * Umlauts and other character modifiers are not given.
        * There is no space after any initial and its following period.
        * Blank spaces are used in a name only if properly part of the surname (e.g., J.VAN DORN),
        or between surname and Jr., II, or III.
        * Abbreviations that are part of a surname, such as Jr., St. or Ste., are followed by a period
        and a space before the next part of the surname.
        * Group names used for one or all of the authors should be spelled out in full.
        * The name of the larger group comes before the name of a subdivision,
        e.g., University of Somewhere, Department of Chemistry.
        * Names are given in English if there is an accepted English version; otherwise in the native language,
        transliterated if necessary.
        """
        return self._authors

    @property
    def pdb_format(self) -> str:
        return

    @property
    def pdb_format_header(self) -> str:
        classification = "/".join(", ".join(entry) for entry in self.classification)
        dep_date = fields.Date.to_pdb(self.deposition_date)
        return f"HEADER{'':4}{classification:<40}{dep_date}{'':3}{self.pdb_id}{'':14}"

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
    Primary structure section of the PDB file, with information on the primary structure of the polymers
    present in the PDB file, including the sequence of residues in each chain of the macromolecule(s),
    and the relationship between the PDB sequence and that found in the sequence database(s).

    Notes
    -----
    * This section contains all data from the PDB records DBREF, DBREF1, DBREF2, SEQADV, SEQRES, and MODRES.
    * Embedded in these records are chain identifiers and sequence numbers that allow other records
    to link into the sequence.
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
        DBREF and DBREF1/DBREF2 records of the PDB file, providing cross-references between
        the sequences of the polymers (chains) in the PDB file (as it appears in the SEQRES
        record), and a corresponding database sequence.

        Returns
        -------
        pandas.DataFrame | None
            Columns:
            * chain_id: Chain identifier of the polymer in the PDB file.
            * residue_num_begin: Initial residue sequence number of the polymer in the PDB file.
            * residue_icode_begin: Initial residue insertion code of the polymer in the PDB file.
            * residue_num_end: Ending residue sequence number of the polymer in the PDB file.
            * residue_icode_end: Ending residue insertion code of the polymer in the PDB file.
            * db: Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            * db_chain_accession: Accession code of the polymer in the database.
            * db_chain_id: Reference to 'chain_id' in the database.
            * db_residue_num_begin: Reference to 'residue_num_begin' in the database.
            * db_residue_icode_begin: Reference to 'residue_icode_begin' in the database.
            * db_residue_num_end: Reference to 'residue_num_end' in the database.
            * db_residue_icode_end: Reference to 'residue_icode_end' in the database.

        Notes
        -----
        * PDB entries contain multi-chain molecules with sequences that may be wild type, variant,
        or synthetic. Sequences may also have been modified through site-directed mutagenesis
        experiments (engineered). A number of PDB entries report structures of individual domains
        cleaved from larger molecules. The DBREF records present sequence correlations between
        SEQRES records in the PDB file, and corresponding GenBank (for nucleic acids) or
        UNIPROT/Norine (for proteins) entries. PDB entries containing heteropolymers are linked
        to different sequence database entries.
        * If no reference is found in the sequence databases, then the PDB entry itself is
        given as the reference.
        * All polymers in the entry must be assigned a DBREF (or DBREF1/DBREF2) record.
        * DBREF and DBREF1/DBREF2 records contain the same information. DBREF1/DBREF2 records
        are a two-line format record, used when the accession code or sequence numbering does
        not fit the space allotted in the standard DBREF format in the PDB file.
        """
        return self._dbref

    @property
    def seqadv(self) -> Optional[pd.DataFrame]:
        """
        SEQADV records of the PDB file, identifying the differences between sequence information
        in the SEQRES records of the PDB entry and the sequence database entry given in DBREF.
        No assumption is made as to which database contains the correct data.

        Returns
        -------
        pandas.DataFrame | None
            Columns:
            * residue_name: Name of the conflicting residue in the PDB file.
            * chain_id: Chain identifier of the conflicting residue's parent polymer in the PDB
            file.
            * residue_num: Sequence number of the conflicting residue in the PDB file.
            * residue_icode: Insertion code of the conflicting residue in the PDB file.
            * db: Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            * db_chain_accession: Accession code of the polymer (chain) in the database.
            * db_residue_name: Reference to 'residue_name' in the database.
            * db_residue_num: Reference to 'residue_num' in the database.
            * description: Description of the conflict.
                * Some possible comments are: 'Cloning artifact', 'Expression tag', 'Conflict', 'Engineered',
                'Variant', 'Insertion', 'Deletion', 'Microheterogeneity', 'Chromophore'.
                * If a conflict is not classifiable by these terms, a reference to either a published paper,
                a PDB entry, or a REMARK within the entry is given. The comment 'SEE REMARK 999' is used
                when the comment is too long.

        Notes
        -----
        * In a number of cases, conflicts between the sequences found in PDB entries and in
        sequence database reference entries have been noted. There are several possible reasons
        for these conflicts, including natural variants or engineered sequences (mutants),
        polymorphic sequences, or ambiguous or conflicting experimental results. These
        discrepancies are reported in this record.
        """
        return self._seqadv

    @property
    def seqres(self) -> Optional[pd.DataFrame]:
        """
        SEQRES records of the PDB file, containing information on the sequence of each
        polymer (i.e. chain) in the file, that is, a listing of the consecutive chemical components
        covalently linked in a linear fashion to form a polymer.

        Returns
        -------
        pandas.DataFrame | None
            The columns are defined as follows:
            chain_id: Chain identifier of the polymer in the PDB file.
            residue_count: Number of residues in the polymer (chain).
            residue_names: Name of the residues in the polymer (chain).

        Notes
        -----
        * The components (i.e. residues) of each sequence may be standard or modified amino/nucleic acids,
        or other residues that are linked to the standard backbone in the polymer.
        Components that are linked to side-chains, or sugars and/or bases are not listed here.
        * Ribo- and deoxyribonucleotides are distinguished; ribo residues are identified with the
        residue names A, C, G, U and I, while deoxy residues are identified with the residue names
        DA, DC, DG, DT and DI. Modified nucleotides are marked by separate 3-letter residue codes.
        * Residues in the ATOM records must agree with the corresponding sequence in SEQRES records.
        * Known problems:
          * Polysaccharides are not properly represented.
          * If the starting position of a sequence is unknown, the sequence cannot be described.
          * For cyclic peptides, a random residue must be assigned as the N-terminus.
        """
        return self._seqres

    @property
    def modres(self) -> Optional[pd.DataFrame]:
        """
        MODRES records of the PDB file, providing descriptions of modifications
        (e.g. chemical or post-translational) to protein and nucleic acid residues.
        Included are correlations between residue names given in a PDB entry and standard residues.

        Returns
        -------
        pandas.DataFrame | None
            Columns:
            * residue_name: Name of the modified residue, as used in the PDB file.
            * chain_id: Chain identifier of the modified residue's parent chain in the PDB file.
            * residue_num: Sequence number of the modified residue in the PDB file.
            * residue_icode: Insertion code of the modified residue in the PDB file.
            * residue_name_std: Standard name of the modified residue.
            * description: Description of the modification.

        Notes
        -----
        * Residues modified post-translationally, enzymatically, or by design are described.
        In those cases where the wwPDB has opted to use a non-standard residue name for the
        residue, MODRES also correlates the new name to the precursor standard residue name.
        * D-amino acids are given their own residue name (resName), i.e., DAL for D-alanine.
        The residue name appears in the SEQRES records, and has the associated MODRES, HET, and FORMUL records.
        The coordinates are given as HETATMs within the ATOM records and occur in the correct order within
        the chain. This ordering is an exception to the Order of Records.
        * When a standard residue name is used to describe a modified site, residue_name and residue_name_std
        contain the same value.
        * MODRES is mandatory when modified standard residues exist in the entry, but is not required if
        coordinate records are not provided for the modified residue.
        """
        return self._modres


class Heterogen:
    """
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
    """
    def __init__(
            self,
            het_data: pd.DataFrame,
            het: pd.DataFrame
    ):
        self._het_data = het_data
        self._het = het
        return

    @property
    def het(self) -> pd.DataFrame:
        """
        HET records of the PDB file. Each non-standard group (residue) is assigned a hetID of
        max. 3 alphanumeric characters. The sequence number, chain identifier, insertion code,
        and number of coordinate records are given for each occurrence of the HET group in the entry.

        Returns
        -------
        pandas.DataFrame
            Columns:
            * het_id: Identifier of the non-standard residue; each unique ID represents a unique molecule.
            * chain_id: Chain identifier of the non-standard residue's parent chain in the PDB file.
            * residue_num: Sequence number of the non-standard residue in the PDB file.
            * residue_icode: Insertion code of the non-standard residue in the PDB file.
            * hetatm_count: Number of HETATM records present in the PDB file corresponding to this molecule.
            * description: Description of the non-standard residue.
        """
        return self._het

    @property
    def het_data(self) -> pd.DataFrame:
        """
        HETNAM, HETSYN and FORMUL records of the PDB file, containing the name, synonyms,
        and chemical formulas of each unique non-standard group in the file.

        Returns
        -------
        pandas.DataFrame
            Index:
                * het_id (hetID): Identifier of the non-standard residue;
                each unique ID represents a unique molecule.
            Columns:
            * component_num: The component number of the heterogen group (see Notes for more info).
            * name: Chemical name of the heterogen group.
            * synonyms: Synonyms of the heterogen group.
            * formula: Chemical formula (plus charge) of the heterogen group.
            * count_in_chain: Number of occurrences of the heterogen group within a chain.
            * count_outer_chain: Number of remaining occurrences of the heterogen group. The sum of
            `count_in_chain` and `count_outer_chain` columns equals to the total number of occurrences of
            the group in the file.

        Notes
        -----
        * PDB entries follow IUPAC/IUB naming conventions to describe groups systematically.
        * The special character '~' is used to indicate superscript in a heterogen name.
        For example: N6 will be listed in the HETNAM section as N~6~, with the ~ character
        indicating both the start and end of the superscript in the name, e.g.,
        N-(BENZYLSULFONYL)SERYL-N~1~-{4-[AMINO(IMINO)METHYL]BENZYL}GLYCINAMIDE.
        * The elements of the chemical formula are given in the order following Hill ordering.
        The order of elements depends on whether carbon is present or not. If carbon is present,
        the order should be: C, then H, then the other elements in alphabetical order of their
        symbol. If carbon is not present, the elements are listed purely in alphabetic order of
        their symbol. This is the 'Hill' system used by Chemical Abstracts.
        * In the chemical formula, the number of each atom type present immediately follows its
        chemical symbol without an intervening blank space. There will be no number indicated
        if there is only one atom for a particular atom type.
        * Each set of SEQRES records and each HET group is assigned a component number in an entry.
        These numbers are assigned serially, beginning with 1 for the first set of SEQRES records.
        In addition:
            * If a HET group is presented on a SEQRES record its FORMUL is assigned the component
            number of the chain in which it appears.
            * If the HET group occurs more than once and is not presented on SEQRES records, the
            component number of its first occurrence is used.
        """
        return self._het_data


class SecondaryStructure:
    """
    Secondary structure section of the PDB file, describing helices and sheets found in the
    protein and polypeptide structures.

    """
    def __init__(self, helix: pd.DataFrame, sheet: pd.DataFrame):
        self._helix: pd.DataFrame = helix
        self._sheet: pd.DataFrame = sheet
        return

    @property
    def helix(self) -> pd.DataFrame:
        """
        HELIX records of the PDB file, describing the helices in the molecule.

        Returns
        -------
        pandas.DataFrame | None
            Index:
            * helix_id (helixID): A unique alphanumeric identifier (max. 3 letters) for each helix.
            Columns:
            * class (helixClass): Classification of the helix as follows:
                * 1: right-handed alpha (default)
                * 2: right-handed omega
                * 3: right-handed pi
                * 4: right-handed gamma
                * 5: right-handed 310
                * 6: left-handed alpha
                * 7: left-handed omega
                * 8: left-handed gamma
                * 9: 27 ribbon/helix
                * 10: polyproline
            * length: Number of residues in the helix.
            * description (comment): Description of the helix.
            * residue_name_begin (initResName): Name of the initial residue (i.e. N-terminal) in the helix.
            * chain_id_begin (initChainID): Chain ID of the initial residue.
            * residue_num_begin (initSeqNum): Residue number of the initial residue.
            * residue_icode_begin (initICode): Insertion code of the initial residue.
            * residue_name_end (endResName): Name of the terminal residue (i.e. C-terminal) in the helix.
            * chain_id_end (endChainID): Chain ID of the terminal residue.
            * residue_num_end (endSeqNum): Residue number of the terminal residue.
            * residue_icode_end (endICode): Insertion code of the terminal residue.
        """
        return self._helix

    @property
    def sheet(self) -> pd.DataFrame:
        """
        SHEET records of the PDB file, describing the sheets in the molecule.

        Returns
        -------
        pandas.DataFrame | None
            Index:
            * sheet_id (helixID): A unique alphanumeric identifier (max. 3 letters) for each helix.
            Columns:
            * class (helixClass): Classification of the helix as follows:
                * 1: right-handed alpha (default)
                * 2: right-handed omega
                * 3: right-handed pi
                * 4: right-handed gamma
                * 5: right-handed 310
                * 6: left-handed alpha
                * 7: left-handed omega
                * 8: left-handed gamma
                * 9: 27 ribbon/helix
                * 10: polyproline
            * length: Number of residues in the helix.
            * description (comment): Description of the helix.
            * residue_name_begin (initResName): Name of the initial residue (i.e. N-terminal) in the helix.
            * chain_id_begin (initChainID): Chain ID of the initial residue.
            * residue_num_begin (initSeqNum): Residue number of the initial residue.
            * residue_icode_begin (initICode): Insertion code of the initial residue.
            * residue_name_end (endResName): Name of the terminal residue (i.e. C-terminal) in the helix.
            * chain_id_end (endChainID): Chain ID of the terminal residue.
            * residue_num_end (endSeqNum): Residue number of the terminal residue.
            * residue_icode_end (endICode): Insertion code of the terminal residue.
        """
        return self._sheet


class Feature:
    """
    Miscellaneous features section of the PDB file, describing properties in the molecule such as
    environments surrounding a non-standard residue or the assembly of an active site.
    """
    def __init__(self, site_data: pd.DataFrame, site_residues: pd.DataFrame):
        self._site_data = site_data
        self._site_residues = site_residues
        return

    @property
    def site_date(self) -> pd.DataFrame:
        """
        REMARK 800 and parts of SITE records of the PDB file, describing each site or environment surrounding
        the ligands in the file. Sites may be of catalytic, co-factor, anti-codon, regulatory or other nature.

        Returns
        -------
        pandas.DataFrame
        """
        return self._site_data

    @property
    def site_residues(self) -> pd.DataFrame:
        """

        Returns
        -------

        """
        return self._site_residues


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