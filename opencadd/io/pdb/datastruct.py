"""
Data structures representing a PDB file.
"""

from typing import Sequence, Optional
import datetime

import numpy as np
import pandas as pd

from opencadd.io.pdb import _fields
from opencadd import _general_exceptions, _typing


class RecordHeader:
    """
    HEADER record of the PDB file, containing the entry's PDB ID, classification, and deposition date.
    """

    def __init__(
            self,
            pdb_id: str,
            dep_date: datetime.date,
            classification: Sequence[Sequence[str]],
    ):
        """
        Parameters
        ----------
        pdb_id : str
            PDB identification code (PDB ID) of the entry.
        dep_date : datetime.date
            Date of deposition of the entry at the Protein Data Bank.
        classification : tuple of tuple of str
            Classification of each molecule within the entry.
            Each sub-tuple corresponds to one molecule in the entry, with each element
            of the sub-tuple describing one classification/function of that molecule.
        """
        self._pdb_id = pdb_id
        self._dep_date = dep_date
        self._classification = classification
        return

    @property
    def pdb_id(self) -> str:
        """
        PDB identification code (PDB ID) of the entry.
        """
        return self._pdb_id

    @property
    def dep_date(self) -> datetime.date:
        """
        Date of deposition of the entry at the Protein Data Bank.
        """
        return self._dep_date

    @property
    def classification(self) -> tuple[tuple[str, ...], ...]:
        """
        Classification of each molecule within the entry.

        Returns
        -------
        tuple of tuple of str
            Each sub-tuple corresponds to one molecule in the entry, with each element
            of the sub-tuple describing one classification/function of that molecule.
        """
        return self._classification

    @pdb_id.setter
    def pdb_id(self, value):
        _fields.IDcode.verify(value)
        self._pdb_id = value
        return

    @dep_date.setter
    def dep_date(self, value):
        _fields.Date.verify(value)
        self._dep_date = value
        return

    @classification.setter
    def classification(self, value):
        if not isinstance(value, _typing.ArrayLike):
            raise TypeError()
        for sub_seq in value:
            if not isinstance(sub_seq, _typing.ArrayLike):
                raise TypeError()
            for elem in sub_seq:
                if not isinstance(elem, str):
                    raise TypeError()
        self._classification = value
        return

    @property
    def pdb_format(self) -> str:
        """
        PDB-formatted string representation of the record.
        """
        classification = "/".join(", ".join(entry) for entry in self.classification)
        dep_date = _fields.Date.to_pdb(self.dep_date)
        return f"HEADER{'':4}{classification:<40}{dep_date}{'':3}{self.pdb_id}{'':14}"

    def __repr__(self):
        return f"Header({self.pdb_id}, {self.dep_date}, {self.classification})"

    def __str__(self):
        lines = [
            f"PDB-ID: {self.pdb_id}",
            f"Deposition Date: {self.dep_date}",
            "Classification (per entity):"
        ] + [f"\t{i+1}. {', '.join(entity)}" for i, entity in enumerate(self.classification)]
        return "\n".join(lines)


class RecordObslte:
    """
    OBSLTE records of the PDB file, indicating the date the entry was removed (“obsoleted”) from the
    PDB's full release, and the PDB IDs of the new entries, if any, that have replaced this entry.

    This record only appears in entries that have been removed from public distribution,
    due to major revisions to coordinates that change the structure's geometry or chemical composition,
    such as changes in polymer sequences, or identity of ligands.
    """
    def __init__(
            self,
            pdb_id: str,
            rep_date: datetime.date,
            rep_pdb_id: Sequence[str],
    ):
        """
        Parameters
        ----------
        pdb_id : str
            PDB identification code (PDB ID) of the entry.
        rep_date : datetime.date
            Date of removal of the entry from the PDB's full release.
        rep_pdb_id : sequence of str
            PDB IDs of the new entries that have replaced this entry.
        """
        self._pdb_id = pdb_id
        self._rep_date = rep_date
        self._rep_pdb_id = np.asarray(rep_pdb_id)
        return

    @property
    def pdb_id(self) -> str:
        """
        PDB identification code (PDB ID) of the entry.
        """
        return self._pdb_id

    @property
    def rep_date(self) -> datetime.date:
        """
        The date the entry was removed (“obsoleted”) from the PDB's full release.
        """
        return self._rep_date

    @property
    def rep_pdb_id(self) -> np.ndarray:
        """
        PDB IDs of the new entries that have replaced this entry.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: <U4]
            1D array of PDB IDs as 4-character strings.
        """
        return self._rep_pdb_id

    @pdb_id.setter
    def pdb_id(self, value):
        _fields.IDcode.verify(value)
        self._pdb_id = value
        return

    @rep_date.setter
    def rep_date(self, value):
        _fields.Date.verify(value)
        self._rep_date = value
        return

    @rep_pdb_id.setter
    def rep_pdb_id(self, value):
        _fields.IDcode.verify(value)
        self._rep_pdb_id = np.asarray(value)
        return

    def __repr__(self):
        return f"Obslte({self.pdb_id}, {self._rep_date}, {self.rep_pdb_id})"

    def __str__(self):
        lines = [
                    f"PDB-ID: {self.pdb_id}",
                    f"Replacement Date: {self.rep_date}",
                    f"Replacement PDB IDs: {', '.join(self.rep_pdb_id)}"
        ]
        return "\n".join(lines)

    @property
    def pdb_format(self) -> str:
        """
        PDB-formatted string representation of the record.
        """
        raise NotImplementedError


class RecordCaveat:
    """
    CAVEAT records of the PDB file, containing a free text description of
    errors and unresolved issues in the entry, if any.

    Notes
    -----
    * This record also appears in entries for which the Protein Data Bank is unable to verify the
      transformation of the coordinates back to the crystallographic cell.
      In these cases, the molecular structure may still be correct.
    """

    def __init__(self, pdb_id: str, description: str):
        """
        Parameters
        ----------
        pdb_id : str
            PDB identification code (PDB ID) of the entry.
        description : str
            A free text, describing the errors and unresolved issues in the entry.
        """
        self.pdb_id = pdb_id
        self.description = description
        return

    @property
    def pdb_id(self) -> str:
        """
        PDB identification code (PDB ID) of the entry.
        """
        return self._pdb_id

    @property
    def description(self) -> str:
        """
        A free text string, describing the errors and unresolved issues in the entry.
        """
        return self._description

    @pdb_id.setter
    def pdb_id(self, value):
        _fields.IDcode.verify(value)
        self._pdb_id = value
        return

    @description.setter
    def description(self, value):
        if not isinstance(value, str):
            raise TypeError()
        self._description = value
        return

    @property
    def pdb_format(self) -> str:
        """
        PDB-formatted string representation of the record.
        """
        raise NotImplementedError

    def __repr__(self):
        return f"Caveat({self.pdb_id}, {self.description})"

    def __str__(self):
        lines = [
                    f"PDB-ID: {self.pdb_id}",
                    f"Caveat: {self.description}",
        ]
        return "\n".join(lines)


class RecordSPRSDE:

    def __init__(self, pdb_id: str, sprsde_date: datetime.date, sprsde_pdb_id: np.ndarray):
        self._pdb_id = pdb_id
        self._sprsde_date = sprsde_date
        self._sprsde_pdb_id = sprsde_pdb_id

    def __repr__(self):
        return f"RecordSPRSDE({self.pdb_id}, {self.sprsde_date}, {self.sprsde_pdb_id})"

    def __str__(self):
        lines = [
                    f"PDB-ID: {self.pdb_id}",
                    f"Superseded Date: {self.sprsde_date}"
                    f"Superseded PDB IDs: {self.sprsde_pdb_id}",
        ]
        return "\n".join(lines)


    @property
    def pdb_id(self) -> str:
        """
        PDB identification code (PDB ID) of the entry.
        """
        return self._pdb_id

    @property
    def sprsde_pdb_id(self) -> Optional[np.ndarray]:
        """
        PDB IDs of entries that were made obsolete by this entry. The date of this event is given
        in the `superseded_date` property of this object (i.e. `Title.superseded_date`).

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: <U4] | None
            PDB IDs of the superseded entries, as an array of strings.

        Notes
        -----
        * This property corresponds to the 'sIdCode' fields of the SPRSDE records in the PDB file.
        """
        return self._sprsde_pdb_id

    @property
    def sprsde_date(self) -> Optional[datetime.date]:
        """
        The date this entry superseded (i.e. replaced) other entries. The list of superseded PDB IDs is given
        in the `superseded_pdb_ids` property of this object (i.e. `Title.superseded_pdb_ids`).

        Returns
        -------
        datetime.date | None
            The date this entry superseded other entries.

        Notes
        -----
        * This property corresponds to the 'sprsdeDate' fields of the SPRSDE records in the PDB file.
        """
        return self._sprsde_date


class RecordSite:
    """

    """
    def __init__(self, site_data: pd.DataFrame, site_residues: pd.DataFrame):
        self._site_data = site_data
        self._site_residues = site_residues
        return

    @property
    def site_data(self) -> pd.DataFrame:
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


class RecordJRNL:
    def __init__(
            self,
            author: Optional[np.ndarray] = None,
            title: Optional[str] = None,
            editor: Optional[np.ndarray] = None,
            pub_name: Optional[str] = None,
            vol: Optional[str] = None,
            page: Optional[str] = None,
            year: Optional[int] = None,
            pub: Optional[str] = None,
            issn: Optional[str] = None,
            essn: Optional[str] = None,
            pm_id: Optional[str] = None,
            doi: Optional[str] = None
    ):
        self._author = author
        self._title = title
        self._editor = editor
        self._pub_name = pub_name
        self._vol = vol
        self._page = page
        self._year = year
        self._pub = pub
        self._issn = issn
        self._essn = essn
        self._pm_id = pm_id
        self._doi = doi
        return

    @property
    def author(self):
        return self._author

    @property
    def title(self):
        return self._title

    @property
    def editor(self):
        return self._editor

    @property
    def pub_name(self):
        return self._pub_name

    @property
    def vol(self):
        return self._vol

    @property
    def page(self):
        return self._page

    @property
    def year(self):
        return self._year

    @property
    def pub(self):
        return self._pub

    @property
    def issn(self):
        return self._issn

    @property
    def essn(self):
        return self._essn

    @property
    def pm_id(self):
        return self._pm_id

    @property
    def doi(self):
        return self._doi

    @property
    def url(self):
        return f"https://doi.org/{self.doi}"

    def __repr__(self):
        arguments = []
        for prop, name in (
                (self.author, "author"),
                (self.title, "title"),
                (self.editor, "editor"),
                (self.pub_name, "pub_name"),
                (self.vol, "vol"),
                (self.page, "page"),
                (self.year, "year"),
                (self.pub, "pub"),
                (self.issn, "issn"),
                (self.essn, "issn"),
                (self.pm_id, "pm_id"),
                (self.doi, "doi")
        ):
            if prop is not None:
                if isinstance(prop, str):
                    arguments.append(f"{name}='{prop}'")
                else:
                    arguments.append(f"{name}={prop}")
        return f"RecordJRNL({', '.join(arguments)})"


class RecordREMARK:
    def __init__(
            self,
            full_text: dict,
            resolution: Optional[float] = None,
            version: Optional[dict] = None,
    ):
        self._full_text = full_text
        self._resolution = resolution
        self._version = version
        return

    def __call__(self, remark_num: int, printout: bool = True):
        if remark_num not in self._full_text:
            return
        remark_lines = self._full_text[remark_num]
        if printout:
            print("\n".join(remark_lines))
        else:
            return remark_lines

    @property
    def resolution(self):
        return self._resolution

    @property
    def version(self):
        return self._version


class RecordSITE:
    def __init__(self, data: pd.DataFrame, residues: pd.DataFrame):
        self._data = data
        self._residues = residues
        return

    @property
    def data(self):
        return self._data

    @property
    def residues(self):
        return self._residues


class RecordCRYST1:
    """
    CRYST1 record of the PDB file, containing the unit cell parameters, space group and Z value.

    Notes
    -----
    * If the structure was not determined by crystallographic means, CRYST1 contains the unitary values, i.e.:
        * a = b = c = 1.0
        * α = β = γ = 90 degrees
        * space group = P 1
        * Z = 1
    """
    def __init__(
            self,
            lengths: np.ndarray,
            angles: np.ndarray,
            z: int,
            space_group: str,
    ):
        self._lengths = lengths
        self._angles = angles
        self._z = z
        self._space_group = space_group
        return

    @property
    def lengths(self) -> np.ndarray:
        """
        Lattice parameters a, b, and c, i.e. the lengths of the unit cell, in Ångstrom.

        Returns
        -------
        numpy.ndarray, shape: (3,), dtype: float64
            array(a, b, c)
        """
        return self._lengths

    @property
    def angles(self) -> np.ndarray:
        """
        Lattice parameters α, β, and γ, i.e. the angles between the edges of the unit cell, in degrees.

        Returns
        -------
        numpy.ndarray, shape: (3,), dtype: float64
            array(α, β, γ)
        """
        return self._angles

    @property
    def z(self) -> int:
        """
        Z-value of the unit cell, i.e. the number of polymeric chains in a unit cell.
        In the case of heteropolymers, Z is the number of occurrences of the most populous chain.

        Returns
        -------
        int

        Notes
        -----
        As an example, given two chains A and B, each with a different sequence, and the space group 'P 2'
        that has two equipoints in the standard unit cell, the following table gives the correct Z value:
            Asymmetric Unit Content     Z value
            -----------------------     -------
                                  A     2
                                 AA     4
                                 AB     2
                                AAB     4
                               AABB     4
        """
        return self._z

    @property
    def space_group(self) -> str:
        """
        Hermnn-Mauguin space group symbol.

        Returns
        -------
        str
            Full international Table's Hermann-Mauguin symbol, without parenthesis, and the screw axis
            described as a two-digit number.
            Examples: 'P 1 21 1' (instead of 'P 21'), 'P 43 21 2'.

        Notes
        -----
        * For a rhombohedral space group in the hexagonal setting, the lattice type symbol used is H.
        """
        return self._space_group

    def __repr__(self):
        return f"RecordCRYST1({self.lengths}, {self.angles}, {self.z}, {self.space_group})"

    def __str__(self):
        return (
            f"Space group: {self.space_group}\n"
            f"Z-value: {self.z}\n"
            f"Unit cell parameters:\n"
            f"\tLengths (a, b, c): {self.lengths}\n"
            f"\tAngles (α, β, γ): {self.angles}"
        )


class RecordXForm:
    def __init__(self, matrix: np.ndarray, vector: np.ndarray):
        self._matrix = matrix
        self._vector = vector
        return

    @property
    def matrix(self):
        return self._matrix

    @property
    def vector(self):
        return self._vector

    def __repr__(self):
        return f"RecordXForm({self.matrix}, {self.vector})"

    def __str__(self):
        return (
            f"Transformation Matrix: {self.matrix}\n"
            f"Translation Vector: {self.vector}"
        )


class RecordMTRIX:
    def __init__(self, serial: np.ndarray, matrices: np.ndarray, vectors: np.ndarray, is_given: np.ndarray):
        self._serial = serial
        self._matrices = matrices
        self._vectors = vectors
        self._is_given = is_given
        self._df = pd.DataFrame({"serial": serial, "is_given": is_given}).set_index("serial", drop=False)
        self._xforms = [RecordXForm(matrix=matrix, vector=vector) for matrix, vector in zip(matrices, vectors)]
        return
    


class PDBFile:
    def __init__(
            self,
            header: Optional[RecordHeader] = None,
            obslte: Optional[RecordObslte] = None,
            title: Optional[str] = None,
            split: Optional[np.ndarray] = None,
            caveat: Optional[str] = None,
            compnd: Optional[pd.DataFrame] = None,
            source: Optional[pd.DataFrame] = None,
            keywds: Optional[np.ndarray] = None,
            expdta: Optional[np.ndarray] = None,
            nummdl: Optional[int] = None,
            mdltyp: Optional[np.ndarray] = None,
            author: Optional[np.ndarray] = None,
            revdat: Optional[np.ndarray] = None,
            sprsde: Optional[RecordSPRSDE] = None,
            jrnl: Optional[RecordJRNL] = None,
            remark: Optional[RecordREMARK] = None,
            dbref: Optional[pd.DataFrame] = None,
            seqadv: Optional[pd.DataFrame] = None,
            seqres: Optional[pd.DataFrame] = None,
            modres: Optional[pd.DataFrame] = None,
            het: Optional[pd.DataFrame] = None,
            hetnam: Optional[pd.DataFrame] = None,
            helix: Optional[pd.DataFrame] = None,
            sheet: Optional[pd.DataFrame] = None,
            ssbond: Optional[pd.DataFrame] = None,
            link: Optional[pd.DataFrame] = None,
            cispep: Optional[pd.DataFrame] = None,
            site: Optional[pd.DataFrame] = None,
            cryst1: Optional[RecordCRYST1] = None,
            origx: Optional[None] = None,
            scale: Optional[None] = None,
            mtrix: Optional[pd.DataFrame] = None,
            atom: Optional[pd.DataFrame] = None,
            anisou: Optional[pd.DataFrame] = None,
            ter: Optional[pd.DataFrame] = None,
            conect: Optional[pd.DataFrame] = None,
    ):
        self._header = header
        self._obslte = obslte
        self._title = title
        self._split = split
        self._caveat = caveat
        self._compnd = compnd
        self._source = source
        self._keywds = keywds
        self._expdta = expdta
        self._nummdl = nummdl
        self._mdltyp = mdltyp
        self._author = author
        self._revdat = revdat
        self._sprsde = sprsde
        self._jrnl = jrnl
        self._remark = remark
        self._dbref = dbref
        self._seqadv = seqadv
        self._seqres = seqres
        self._modres = modres
        self._het = het
        self._hetnam = hetnam
        self._helix = helix
        self._sheet = sheet
        self._ssbond = ssbond
        self._link = link
        self._cispep = cispep
        self._site = site
        self._cryst1 = cryst1
        self._origx = origx
        self._scale = scale
        self._mtrix = mtrix
        self._atom = atom
        self._anisou = anisou
        self._ter = ter
        self._conect = conect
        return

    @property
    def header(self) -> Optional[RecordHeader]:
        """
        HEADER record of the PDB file, containing the entry's PDB ID, classification, and deposition date.

        Returns
        -------
        Header or None
            If the PDB file contains no HEADER record, `None` is returned,
            otherwise an instance of `opencadd.io.pdb.datastruct.RecordHeader` with following properties:

            pdb_id : str
                PDB identification code (PDB ID) of the entry.
            dep_date : datetime.date
                Date of deposition of the entry at the Protein Data Bank.
            classification : tuple of tuple of str
                Classification of each molecule within the entry.
        """
        return self._header

    @property
    def obslte(self) -> Optional[RecordObslte]:
        """
        OBSLTE records of the PDB file, indicating the date the entry was removed (“obsoleted”) from the
        PDB's full release, and the PDB IDs of the new entries, if any, that have replaced this entry.

        This record only appears in entries that have been removed from public distribution,
        due to major revisions to coordinates that change the structure's geometry or chemical composition,
        such as changes in polymer sequences, or identity of ligands.

        Returns
        -------
        RecordObslte or None
            If the PDB file contains no OBSLTE record, `None` is returned,
            otherwise an instance of `opencadd.io.pdb.datastruct.RecordObslte`.
        """
        return self._obslte

    @property
    def title(self) -> Optional[str]:
        """
        TITLE records of the PDB file, containing a title for the experiment or analysis that is represented
        in the entry.

        The title is a free text, describing the contents of the entry and any procedures or
        conditions that distinguish it from similar entries.
        Some data that may be included are experiment type, description of the mutation,
        and the fact that only alpha carbon coordinates have been provided in the entry.

        Returns
        -------
        str or None
            If the PDB file contains no TITLE record, `None` is returned,
            otherwise the title as a free-text string.
        """
        return self._title

    @property
    def split(self) -> Optional[np.ndarray]:
        """
        SPLIT records of the PDB file, containing the PDB IDs of entries that are required
        to reconstitute a complete complex.

        This record only appears in entries that compose a part of a larger macromolecular complex.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: <U4] or None
            If the PDB file contains no SPLIT record, `None` is returned,
            otherwise a 1D array of PDB IDs as 4-character strings.
        """
        return self._split

    @property
    def caveat(self) -> Optional[RecordCaveat]:
        """
        CAVEAT records of the PDB file, containing a free text description of
        errors and unresolved issues in the entry, if any.

        Returns
        -------
        Caveat or None
            If the PDB file contains no CAVEAT record, `None` is returned,
            otherwise an instance of `opencadd.io.pdb.datastruct.record.Caveat`.

        Notes
        -----
        * This record also appears in entries for which the Protein Data Bank
        """
        return self._caveat

    @property
    def compnd(self) -> Optional[pd.DataFrame]:
        """
        COMPND records of the PDB file, describing the macromolecular contents of the PDB file,
        or a standalone drug or inhibitor in cases where the entry does not contain a polymer.

        Returns
        -------
        pandas.DataFrame or None
            If the PDB file contains no COMPND record, `None` is returned,
            otherwise a `DataFrame` with columns:

            mol_id (index) : int
                Enumerates each molecule; the same ID appears also in the SOURCE records.
            name : str
                Name of the (macro)molecule. For chimeric proteins, the protein name is
                comma-separated and may refer to the presence of a linker, e.g. "protein_1, linker, protein_2".
            chain_ids : numpy.ndarray[dtype: <U1]
                Chain identifiers in the macromolecule.
            fragment : str
                Name or description of a domain or region of the molecule.
            synonyms : numpy.ndarray[dtype: str]:
                Synonyms for the molecule's name.
            enzyme_commission_num : numpy.ndarray[dtype: str]
                Enzyme commision (EC) numbers associated with the molecule.
            engineered : bool
                Whether the molecule was produced using recombinant technology or by purely chemical synthesis.
            mutation : bool
                Whether there is a mutation in the molecule.
            description : str
                Additional free-text comment.

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
        return self._compnd

    @property
    def source(self) -> Optional[pd.DataFrame]:
        """
        SOURCE records of the PDB file, containing information on the biological/chemical source of
        each biological molecule in the PDB file, or a standalone drug or inhibitor in cases
        where the entry does not contain a polymer.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no SOURCE record, `None` is returned,
            otherwise a dataframe with columns:

            mol_id (index) : int
                Enumerates each molecule; the same ID appears also in the `compound` property of
                this object (i.e. `Title.compound`).
            synthetic : str
                Indicates a chemically synthesized source.
            fragment : str
                Specifies a domain or fragment of the molecule.
            organism : str
                Common name of the organism
            organism_sci : str
                Scientific name of the organism.
            organism_tax_id : str
                NCBI Taxonomy ID of the organism.
            strain : str
                Identifies the strain.
            variant : str
                Identifies the variant.
            cell_line : str
                The specific line of cells used in the experiment.
            atcc_id : str
                American Type Culture Collection tissue culture number.
            organ : str
                Organized group of tissues that carries on a specialized function.
            tissue : str
                Organized group of cells with a common function and structure.
            cell : str
                Identifies the particular cell type.
            organelle : str
                Organized structure within a cell.
            secretion : str
                Identifies the secretion, such as saliva, urine, or venom, from which the molecule
                was isolated.
            cell_loc : str
                Identifies the location inside/outside the cell, where the compound was found.
                Examples are: 'extracellular', 'periplasmic', 'cytosol'.
            plasmid : str
                Identifies the plasmid containing the gene.
            gene : str
                Identifies the gene.
            expsys : str
                Expression system, i.e. common name of the organism in which the molecule was expressed.
            expsys_sci : str
                Scientific name of the expression system.
            expsys_tax_id : str
                NCBI Taxonomy ID of the expression system.
            expsys_strain : str
                Strain of the organism in which the molecule was expressed.
            expsys_variant : str
                Variant of the organism used as the expression system.
            expsys_cell_line : str
                The specific line of cells used as the expression system.
            expsys_atcc_id : str
                American Type Culture Collection tissue culture number of the expression system.
            expsys_organ : str
                Specific organ which expressed the molecule.
            expsys_tissue : str
                Specific tissue which expressed the molecule.
            expsys_cell : str
                Specific cell type which expressed the molecule.
            expsys_organelle : str
                Specific organelle which expressed the molecule.
            expsys_cell_loc : str
                Identifies the location inside or outside the cell which expressed the molecule.
            expsys_vector_type : str
                Identifies the type of vector used, i.e. plasmid, virus, or cosmid.
            expsys_vector : str
                Identifies the vector used.
            expsys_plasmid : str
                Plasmid used in the recombinant experiment.
            expsys_gene : str
                Name of the gene used in recombinant experiment.
            details : str
                Other details about the source.

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
    def keywds(self) -> Optional[np.ndarray]:
        """
        KEYWDS records of the PDB file, containing keywords/terms relevant to the PDB file,
        similar to that found in journal articles.
        The provided terms may for example describe functional classification, metabolic role, known
        biological or chemical activity, or structural classification.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str]
            If the PDB file contains no KEYWDS record, `None` is returned, otherwise the keywords
            as an array of strings.

        Notes
        -----
        * The classifications given in `PDBFile.header.classification` are also repeated here,
          with two differences: Unlike in `classification`, here the keywords are not grouped per molecule,
          but they are given unabbreviated.
        * This property corresponds to the 'keywds' fields of the KEYWDS records in the PDB file.
        """
        return self._keywds

    @property
    def expdta(self) -> Optional[np.ndarray]:
        """
        EXPDTA records of the PDB file, identifying the experimental technique used for determining the
        structure. This may refer to the type of radiation and sample, or include the spectroscopic or modeling
        technique.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str] | None
            Array of strings, containing one or several of following allowed values:
            'X-RAY DIFFRACTION', 'FIBER DIFFRACTION', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY',
            'ELECTRON MICROSCOPY', 'SOLID-STATE NMR', 'SOLUTION NMR', 'SOLUTION SCATTERING'.

        Notes
        -----
        * Since October 15, 2006, theoretical models are no longer accepted for deposition. Any
          theoretical models deposited prior to this date are archived at:
          ftp://ftp.wwpdb.org/pub/pdb/data/structures/models
        * This property corresponds to the 'technique' fields of the EXPDATA records in the PDB file.
        """
        return self._expdta

    @property
    def nummdl(self) -> Optional[int]:
        """
        NUMMDL record of the PDB file, indicating the total number of models in the entry.

        Returns
        -------
        int | None
            Total number of models in the PDB file.

        Notes
        -----
        * This property corresponds to the 'modelNumber' field of the NUMMDL record in the PDB file.
        """
        return self._nummdl

    @property
    def mdltyp(self) -> Optional[np.ndarray]:
        """
        MDLTYP records of the PDB file, containing additional structural annotations on the coordinates
        in the PDB file, used to highlight certain features.

        Additional structural annotations pertinent to the coordinates in the PDB file, used to highlight
        certain features.

        Returns
        -------
        numpy.ndarray[ndim: 1, dtype: str] | None
            Array of strings corresponding to a list of annotations.

        Notes
        -----
        *  For entries that are determined by NMR methods and the coordinates deposited are either a
          minimized average or regularized mean structure, the tag "MINIMIZED AVERAGE" will be present as the
          first element of the returned array.
        * Where the entry contains entire polymer chains that have only either C-alpha (for proteins) or
          P atoms (for nucleotides), the contents of such chains will be described along with the
          chain identifier, e.g. " CA ATOMS ONLY, CHAIN A, B". For these polymeric chains,
          REMARK 470 (Missing Atoms) will be omitted.
        * This property corresponds to the 'comment' fields of the MDLTYP record in the PDB file.
        """
        return self._mdltyp

    @property
    def author(self) -> Optional[np.ndarray]:
        """
        AUTHOR records of the PDB file, containing the names of the persons responsible for the contents
        of the entry.

        Returns
        -------
        numpy.ndarray[ndim = 1, dtype = str] | None
            Array of strings corresponding to a list of authors.

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
        * This property corresponds to the 'authorList' fields of the AUTHOR record in the PDB file.
        """
        return self._author

    @property
    def revdat(self) -> Optional[pd.DataFrame]:
        """
        REVDAT records of the PDB file, containing a history of modifications made to the entry since its
        release.

        Returns
        -------
        pandas.DataFrame | None
            mod_num (index) : int
                Enumerates each release/modification, starting at 1 for the initial release.
            date : datetime.date
                Date of release/modification.
            pdb_id : str
                PDB ID of the entry for the specific modification/release.
            is_initial : bool
                Indicating the initial release of the entry. The value is `True` for the row at
                index (mod_num) 1, and `False` for all other rows.
            details : numpy.ndarray[ndim: 1, dtype: str]
                Details of the modification as an array of keywords, which are typically PDB record names
                such as 'JRNL', 'SOURCE', 'TITLE', 'COMPND' etc. The keyword 'VERSN' indicates that the file
                has undergone a change in version; The current version is specified in REMARK 4.

        Notes
        -----
        * This property corresponds to the entire REVDAT record in the PDB file.
        """
        return self._revdat

    @property
    def sprsde(self):
        """
        SPRSDE records of the PDB file, containing a list of PDB IDs of the entries that were made obsolete
        by this entry, and the corresponding dates.

        Returns
        -------
        dict[str, str | datetime.date | numpy.ndarray[ndim: 1, dtype: <U4]] | None
            If the PDB file contains no SPRSDE record, `None` is returned, otherwise a dictionary with keys:
            pdb_id : str
                PDB ID of this entry.
            date : datetime.date
                The date this entry superseded the listed entries.
            superseded_pdb_ids : numpy.ndarray[ndim: 1, dtype: <U4]
                PDB IDs of the superseded entries, as an array of strings.
        """
        return self._sprsde

    @property
    def jrnl(self):
        return self._jrnl

    @property
    def remark(self):
        return self._remark

    @property
    def dbref(self) -> Optional[pd.DataFrame]:
        """
        DBREF and DBREF1/DBREF2 records of the PDB file, providing cross-references between each
        sequence (chain) of the polymers in the PDB file (as it appears in the SEQRES records),
        and corresponding GenBank (for nucleic acids) or UNIPROT/Norine (for proteins) database
        sequence entries.

        PDB entries containing heteropolymers are linked to different sequence database entries.
        If no reference is found in the sequence databases, then the PDB entry itself is given as
        the reference.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains neither DBREF nor DBREF1/DBREF2 records, `None` is returned,
            otherwise a dataframe with columns:

            chain_id (index) : str
                Chain identifier of the polymer in the PDB file.
            residue_num_begin : int
                Initial residue sequence number of the polymer in the PDB file.
            residue_icode_begin : str
                Initial residue insertion code of the polymer in the PDB file.
            residue_num_end : int
                Ending residue sequence number of the polymer in the PDB file.
            residue_icode_end : str
                Ending residue insertion code of the polymer in the PDB file.
            db : str
                Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            db_chain_accession : str
                Accession code of the polymer in the database.
            db_chain_id : str
                Reference to 'chain_id' in the database.
            db_residue_num_begin : int
                Reference to 'residue_num_begin' in the database.
            db_residue_icode_begin : str
                Reference to 'residue_icode_begin' in the database.
            db_residue_num_end : int
                Reference to 'residue_num_end' in the database.
            db_residue_icode_end : str
                Reference to 'residue_icode_end' in the database.

        Notes
        -----
        * PDB entries contain multi-chain molecules with sequences that may be wild type, variant,
          or synthetic. Sequences may also have been modified through site-directed mutagenesis
          experiments (engineered). A number of PDB entries report structures of individual domains
          cleaved from larger molecules.
        * This property corresponds to the DBREF and DBREF1/DBREF2 records in the PDB file, which contain
          the same type of information; DBREF1/DBREF2 records are a two-line format record, used when
          the accession code or sequence numbering does not fit the space allotted in the standard DBREF format.
        * All polymers in the entry must be assigned a database reference.

        * Both DBREF and DBREF1/DBREF2 records contain the same type of information; DBREF1/DBREF2 records
          are a two-line format record, used when the accession code or sequence numbering does not fit
          the space allotted in the standard DBREF format.
        """
        return self._dbref

    @property
    def seqadv(self) -> Optional[pd.DataFrame]:
        """
        SEQADV records of the PDB file, identifying the differences between sequence information
        in the SEQRES records of the PDB entry and the sequence database entry given in DBREF.
        No assumption is made as to which database contains the correct data.

        In a number of cases, conflicts between the sequences found in PDB entries and in
        sequence database reference entries have been noted. There are several possible reasons
        for these conflicts, including natural variants or engineered sequences (mutants),
        polymorphic sequences, or ambiguous or conflicting experimental results. These
        discrepancies are reported in this record.

        Returns
        -------
        pandas.DataFrame | None
            If the PDB file contains no SEQADV record, `None` is returned, otherwise a dataframe with columns:

            chain_id (index) : str
                Chain identifier of the conflicting residue's parent polymer in the PDB file.
            pdb_id : str
                PDB ID of the entry.
            residue_name : str
                Name of the conflicting residue in the PDB file.
            residue_num : int
                Sequence number of the conflicting residue in the PDB file.
            residue_icode : str
                Insertion code of the conflicting residue in the PDB file.
            db : str
                Database name (GB (GenBank), PDB (Protein Data Bank), UNP (UNIPROT), NORINE, UNIMES)
            db_chain_accession : str
                Accession code of the polymer (chain) in the database.
            db_residue_name : str
                Reference to 'residue_name' in the database.
            db_residue_num : int
                Reference to 'residue_num' in the database.
            description : str
                Description of the conflict. Some possible comments are:
                'Cloning artifact', 'Expression tag', 'Conflict', 'Engineered', 'Variant',
                'Insertion', 'Deletion', 'Microheterogeneity', 'Chromophore'.
                If a conflict is not classifiable by these terms, a reference to either a published paper,
                a PDB entry, or a REMARK within the entry is given. The comment 'SEE REMARK 999' is used
                when the comment is too long.

        Notes
        -----
        * Microheterogeneity is to be represented as a variant with one of the possible residues in the site
          being selected (arbitrarily) as the primary residue. The residues that do not match the UNP
          reference are listed with the description 'Microheterogeneity'.
        * This property corresponds to the SEQADV records in the PDB file.
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
    def hetnam(self) -> pd.DataFrame:
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
        return self._hetnam

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

    @property
    def ssbond(self):
        return self._ssbond

    @property
    def link(self):
        return self._link

    @property
    def cispep(self):
        return self._cispep

    @property
    def site(self):
        return self._site

    @property
    def cryst1(self):
        return self._cryst1

    @property
    def origx(self):
        return self._origx

    @property
    def scale(self):
        return self._scale

    @property
    def mtrix(self):
        return self._mtrix

    @property
    def atom(self):
        return self._atom

    @property
    def anisou(self):
        return self._anisou

    @property
    def ter(self):
        return self._ter

    @property
    def conect(self):
        return self._conect

    def to_pdb(self):
        def record_end() -> str:
            """
            END record of a PDB file.
            The END record marks the end of the PDB file, and must appear as the final record in every
            file.
            """
            return f"{'END':<80}"

    @property
    def pdb_format_title(self) -> str:
        """
        TITLE Record formatted as in a PDB file.
        """
        num_lines_needed = int(np.ceil(len(self.title) / 70))
        continuation = _fields.Continuation.to_pdb(num_lines=num_lines_needed)
        title = [self.title[i:i + 70] for i in range(0, len(self._title), 70)]
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