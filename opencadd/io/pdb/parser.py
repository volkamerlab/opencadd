import datetime
import re
from typing import Literal, NoReturn, Optional, Dict, Any, Union, Tuple
import warnings

import numpy as np
import pandas as pd

from opencadd.io import parsing
from opencadd.io._helpers_sysio import filelike_to_data_string
from opencadd.io.pdb import sections, records, fields
from opencadd.api.web import rcsb


class PDBParsingError(Exception):
    pass


class PDBParser:
    """
    rec_names : ndarray, shape: (n,), dtype: <U6
        1-dimensional array of strings, containing the record name of each line in the PDB file,
        stripped of whitespaces.
        For example: array(['HEADER', 'TITLE', 'TITLE', ..., 'CONECT', 'MASTER', 'END'])
    rec_chars : ndarray, shape: (n, m), dtype: <U1
        2-dimensional array of characters, where the first axis is the rows (lines) and second
        axis the columns (characters in a line) in a PDB file. In a standard PDB file,
        the number of columns 'm' must be 80.

    """
    def __init__(self, file, strictness: bool = True):
        """
        Parameters
        ----------
        file : FileLike
            PDB file to parse.
        strictness : Literal[0, 1, 2, 3], optional, default: 0
            Level of strictness for raising exceptions and warnings when encountering mistakes
            in the PDB file:
                0: raise only fatal errors and don't show any warnings.
                1: raise only fatal errors. All other errors are reported as warnings.
                2: raise fatal errors and mistakes resulting in ambiguous data. Inconsequential mistakes are
                reported as warnings.
                3: completely validate the PDB file and raise all errors.
        """
        # Declare attributes
        self._strictness: Literal[0, 1, 2, 3]
        self._lines: np.ndarray
        self._lines_chars: np.ndarray
        """Character view of lines, i.e. 2-d array where first axis is lines and second axis is 
        characters"""
        self._lines_rec_names: np.ndarray
        """Extract all record types as 1-d array of stripped record-name strings."""
        self._lines_rec_content: np.ndarray
        self._mask_record: dict
        self._idx_record: dict
        self._has_record: dict
        self._count_record: dict

        # Read file and set attributes
        self.strictness = strictness
        self._lines = np.array(filelike_to_data_string(file=file).splitlines())
        self._lines_chars = self._lines.view(dtype=(str, 1)).reshape(self._lines.size, -1)
        self._lines_rec_names = np.char.strip(self._lines_chars[:, :6].view(dtype=(str, 6)).reshape(-1))
        self._lines_rec_content = self._lines_chars[:, 6:].view(dtype=(str, 74)).reshape(-1)
        self._mask_record = dict()
        self._idx_record = dict()
        self._has_record = dict()
        self._count_record = dict()
        return

    @property
    def file(self) -> PDBFile:
        pdb_file = PDBFile(
            title=self.section_title,
            remark=self.section_remark,
            primary_structure=self.section_primary_structure,
            heterogen=self.section_heterogen,
            secondary_structure=self.section_secondary_structure,
            connectivity=self.section_connectivity,
            feature=self.section_feature,
            crystallographic=self.section_crystallographic,
            coordinate=self.section_coordinate,
            bookkeeping=self.section_bookkeeping,
        )
        return pdb_file

    @property
    def section_title(self) -> sections.Title:
        header = self.header
        obslte = self.obslte
        if header is None:
            header = (None, None, None)
        if obslte is None:
            obslte = (None, None)
        title = sections.Title(
            pdb_id=header[0],
            deposition_date=header[1],
            classification=header[2],
            replacement_date=obslte[0],
            replacement_pdb_ids=obslte[1],
            title=self.title,
            split_pdb_ids=self.split,
            caveat=self.caveat,
            compound=self.compnd,
            source=self.source,
            keywords=self.keywds,
            experimental_techniques=self.expdta,
            count_models=self.nummdl,
            structural_annotations=self.mdltyp,
            authors=self.author,
        )
        return title

    @property
    def section_remark(self) -> sections.Remark:
        pass

    @property
    def section_primary_structure(self) -> sections.PrimaryStructure:
        primary_structure = sections.PrimaryStructure(
            dbref=self.dbref,
            seqadv=self.seqadv,
            seqres=self.seqres,
            modres=self.modres
        )
        return primary_structure

    @property
    def section_heterogen(self) -> sections.Heterogen:
        het_data = pd.concat((self.hetnam, self.hetsyn, self.formul), axis=1)
        heterogen = sections.Heterogen(
            het_data=het_data[
                ["component_num", "name", "is_water", "formula", "count_in_chain", "count_outer_chain"]
            ],
            het=self.het
        )
        return heterogen

    def section_secondary_structure(self) -> sections.SecondaryStructure:
        secondary_structure = sections
        pass

    def section_connectivity_annotation(self) -> sections.Connectivity:
        return

    def section_feature(self):
        pass

    def section_crystallographic(self) -> sections.Crystallographic:
        pass

    def section_coordinate_transformation(self):
        pass

    def section_coordinate(self):
        if not (self.has_record(Record.ATOM) or self.has_record(Record.HETATM)):
            return None
        df_ter_records = self.ter
        mask_lines = np.logical_or(
            self.mask_record(Record.ATOM), self.mask_record(Record.HETATM)
        )
        if not self.has_record(Record.MODEL):
            lines = self._lines_chars[mask_lines]
            atom_data = parsing.extract_columns(lines, Record.ATOM.columns)
            atom_data["model_nr"] = np.full(shape=lines.shape[0], fill_value=1, dtype=np.ubyte)
            atom_data["is_standard_residue"] = self._lines_rec_names[mask_lines] == "ATOM"
            atom_data["is_polymer"] = atom_data["serial"] < df_ter_records.iloc[-1]["serial"]
            return pd.DataFrame(atom_data)

        num_coordinates = np.count_nonzero(mask_lines)

        if self.has_record(Record.MODEL):
            indices_model_start = self.indices_record_lines(Record.MODEL)
            indices_model_end = self.indices_record_lines(Record.ENDMDL)
            for idx_start, idx_end in zip(indices_model_start, indices_model_end):
                pass


        # df_anisou = parsing.dataframe_from_string_array(
        #     self._rec_chars[self._mask_record[Record.ANISOU]],
        #     column_data=Record.ANISOU.columns
        # )

        df_atom_data = parsing.dataframe_from_string_array(
            self._lines_chars[mask_lines],
            column_data=Record.ATOM.columns
        )
        is_polymer = np.zeros(shape=df_atom_data.index.values.size, dtype=np.bool_)
        if self.has_record(Record.MODEL):
            size_models = np.array(
                [
                    np.count_nonzero(mask_lines[start:end])
                    for start, end in zip(
                        self.indices_record_lines(Record.MODEL), self.indices_record_lines(Record.ENDMDL)
                    )
                ]
            )
            df_atom_data["model_nr"] = np.repeat(
                np.arange(1, size_models.size + 1), repeats=size_models
            )
            if self.has_record(Record.TER):
                idx_ter = self.indices_record_lines(Record.TER)
                for idx_model in range(size_models.size):
                    model_start_line = self.indices_record_lines(Record.MODEL)[idx_model]
                    model_end_line = self.indices_record_lines(Record.ENDMDL)[idx_model]
                    ter_mask = np.logical_and(
                        idx_ter > model_start_line,
                        idx_ter < model_end_line
                    )
                    num_polymer_atoms = np.count_nonzero(
                        mask_lines[model_start_line:idx_ter[ter_mask][-1]]
                    )
                    start_idx = np.sum(size_models[:idx_model])
                    is_polymer[start_idx:start_idx+num_polymer_atoms] = True
        else:
            df_atom_data["model_nr"] = 1
            if self.has_record(Record.TER):
                idx_last_ter = self.indices_record_lines(Record.TER)[-1]
                is_polymer[: -np.count_nonzero(mask_lines[idx_last_ter + 1:])] = True
        df_atom_data["is_polymer"] = is_polymer
        return df_atom_data

    def section_connectivity(self):
        pass

    def section_bookkeeping(self):
        pass

    @property
    def header(self) -> Optional[Tuple[str, datetime.date, Tuple[Tuple[str]]]]:
        """
        HEADER record of the PDB file, containing the PDB ID, classification, and deposition date of the entry.

        Returns
        -------
        tuple[str, datetime.date, tuple[tuple[str]]] | None
            PDB ID, deposition date, and classification.
            For the classification, each sub-tuple corresponds to one molecule in the entry, with each element
            of the sub-tuple describing one classification/function of the molecule.
            If the PDB file contains no HEADER record, `None` is returned.

        Notes
        -----
        * Classification may be based on function, metabolic role, molecule type, cellular location
        etc., but it exactly matches one of a collection of strings, available at:
        https://www.wwpdb.org/docs.html
        *  Due to the limited length of the classification field, strings must sometimes be
        abbreviated. In these cases, the full terms are given in KEYWDS records in no strict order.
        """
        if not self.has_record(Record.HEADER):
            if self.strictness:
                self.header_validation()
            return
        # HEADER is a one-time/single-line record; take the first occurrence:
        data = self.extract_record_fields(Record.HEADER)
        if self.strictness:
            self.header_validation(data)
        # The classification string is left-justified, and can describe dual functions of
        # molecules (when applicable) separated by a comma “,”. Entries with multiple
        # molecules in a complex will list the classifications of each macromolecule
        # separated by slash “/”.
        # First, split by '/' and then by ',' to get a tuple of tuples:
        classification = data["classification"][0]
        class_per_entity = tuple(classification.split("/"))
        class_per_entity_and_function = tuple(
            tuple(entity_function.strip() for entity_function in entity_class.split(","))
            for entity_class in class_per_entity
        )
        return data["pdb_id"][0], data["deposition_date"][0], class_per_entity_and_function

    def header_validation(self, fields: dict = None):
        if not self._validate_existence_mandatory_record(Record.HEADER):
            return False
        validated = self._validate_one_time_single_line_record(Record.HEADER)
        validated = validated and self.validate_empty_columns(Record.HEADER)
        if fields is not None:
            validated = validated and self._validate_pdb_ids(fields["pdb_id"])

        # Check position in file; HEADER must always be in first line:
        indices_lines = self.indices_record_lines(Record.HEADER)
        if indices_lines[0] != 0:
            validated = False
            self._raise_or_warn(
                f"HEADER record found at line {indices_lines[0]}. "
                f"It must always be the first line of a PDB file.",
                raise_level=3
            )
        return validated

    @property
    def obslte(self) -> Optional[Tuple[datetime.date, np.ndarray]]:
        """
        OBSLTE record of the PDB file, indicating the date the entry was removed (“obsoleted”) from the
        PDB's full release, and the PDB IDs of the new entries, if any, that have replaced this entry.
        This record only appears in entries that have been removed from public distribution.

        Returns
        -------
        tuple[datetime.date, numpy.ndarray[ndims:1, dtype:<U4]] | None
            Replacement date, and the PDB-IDs that have replaced this entry.
            If the PDB file contains no OBSLTE record, `None` is returned.

        Notes
        -----
        * All OBSLTE entries are available from the PDB archive at:
        https://ftp.wwpdb.org/pub/pdb/data/structures/obsolete.
        """
        if not self.has_record(Record.OBSLTE):
            return
        # OBSLTE is a one-time/multiple-lines record
        fields = self.extract_record_fields(Record.OBSLTE)
        replaced_pdb_ids = fields["replaced_pdb_ids"][fields["replaced_pdb_ids"] != ""]
        return fields["replacement_date"][0], replaced_pdb_ids

    @property
    def title(self) -> Optional[str]:
        """
        TITLE record of the PDB file, containing a title for the experiment or analysis that is represented
        in the entry. Some data that may be included are e.g. experiment type, description of the mutation,
        and the fact that only alpha carbon coordinates have been provided in the entry.

        Returns
        -------
        str | None
            A free text describing the contents of the entry and any procedures or conditions that
            distinguish it from similar entries.
            If the PDB file contains no TITLE record, `None` is returned.
        """
        if not self.has_record(Record.TITLE):
            return
        data = self.extract_record_fields(Record.TITLE)
        if self.strictness:
            self.title_validation(data)
        # TITLE is a one-time/multiple-lines record
        return data["title"]

    def title_validation(self, fields: dict = None):
        pass

    @property
    def split(self) -> Optional[np.ndarray]:
        """
        SPLIT record of the PDB file, containing the PDB IDs of entries that are required to reconstitute a
        complete complex. This only appears in instances where the entry composes a part of a large
        macromolecular complex.

        Returns
        -------
        numpy.ndarray[ndims:1, dtype:<U4] | None
            PDB IDs that are required to reconstitute a complete complex.
            If the PDB file contains no SPLIT record, `None` is returned.

        Notes
        -----
        * If SPLIT records exist, REMARK 350 will contain an amended statement to reflect the entire complex.
        """
        if not self.has_record(Record.SPLIT):
            return
        # SPLIT is a one-time/multiple-lines record
        fields = self.extract_record_fields(Record.SPLIT)
        pdb_ids = fields["pdb_ids"][fields["pdb_ids"] != ""]
        return pdb_ids

    @property
    def caveat(self) -> Optional[str]:
        """
        CAVEAT record of the PDB file, containing a free text comment that warns of errors and unresolved
        issues in the entry. It will also be included in cases where the PDB is unable to verify the
        transformation of the coordinates back to the crystallographic cell. In these cases, the molecular
        structure may still be correct.

        Returns
        -------
        str | None
            A free text comment on the errors and unresolved issues in the entry.
            If the PDB file contains no CAVEAT record, `None` is returned.
        """
        if not self.has_record(Record.CAVEAT):
            return
        # CAVEAT is a one-time/multiple-lines record
        return self.extract_record_fields(Record.CAVEAT)["description"]

    @property
    def compnd(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.COMPND):
            return


        def parse_token_value_pairs(start: int, end: int) -> dict[str, Any]:
            data = dict()
            for token, value in token_value_pairs[start:stop]:
                try:
                    col_name, parser = correspondence[token]
                except KeyError:
                    self._raise_or_warn(
                        f"Encountered invalid token in COMPND record: {token}.",
                        raise_level=2
                    )
                else:
                    data[col_name] = parser(value)
            return data

        df = pd.DataFrame(columns=[col_name for col_name, _ in correspondence.values()])
        token_value_pairs = self.extract_record_fields(Record.COMPND)["compound"]
        ind_mols = np.argwhere(token_value_pairs[:, 0] == "MOL_ID").reshape(-1)
        for start, stop in np.lib.stride_tricks.sliding_window_view((*ind_mols, None), 2):
            mol_id = int(token_value_pairs[start, 1])
            ind_fragments = np.argwhere(token_value_pairs[start+1:stop, 0] == "FRAGMENT").reshape(-1)
            if ind_fragments.size <= 1:
                df.loc[mol_id] = parse_token_value_pairs(start+1, stop)
            else:
                data_mol = parse_token_value_pairs(start+1, start+1+ind_fragments[0])
                for start_frag, end_frag in np.lib.stride_tricks.sliding_window_view((*ind_fragments, None), 2):
                    data_frag = parse_token_value_pairs(start+1+start_frag, start+1+end_frag)
                    sub_df = pd.DataFrame(data_mol | data_frag, index=[mol_id])
                    df = pd.concat((df, sub_df))
        df.index.name = "mol_id"
        return df

    @property
    def source(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.SOURCE):
            return
        df = pd.DataFrame(columns=list(correspondence.values()))
        token_value_pairs = self.extract_record_fields(Record.SOURCE)["source_name"]
        ind_mols = np.argwhere(token_value_pairs[:, 0] == "MOL_ID").reshape(-1)
        index, rows = list(), list()
        for start, stop in np.lib.stride_tricks.sliding_window_view((*ind_mols, None), 2):
            mol_id = int(token_value_pairs[start, 1])
            data = dict()
            for token, value in token_value_pairs[start + 1:stop]:
                try:
                    col_name = correspondence[token]
                except KeyError:
                    self._raise_or_warn(
                        f"Encountered invalid token in SOURCE record: {token},",
                        raise_level=2
                    )
                else:
                    data[col_name] = value
            index.append(mol_id)
            rows.append(data)
        df = pd.DataFrame(rows, index=index)
        df.index.name = "mol_id"
        return df

    @property
    def keywds(self) -> Optional[np.ndarray]:
        if not self.has_record(Record.KEYWDS):
            return
        return self.extract_record_fields(Record.KEYWDS)["keywords"]

    @property
    def expdta(self) -> Optional[np.ndarray]:
        if not self.has_record(Record.EXPDTA):
            return
        return self.extract_record_fields(Record.EXPDTA)["technique"]

    @property
    def nummdl(self) -> Optional[int]:
        if not self.has_record(Record.NUMMDL):
            return
        return self.extract_record_fields(Record.NUMMDL)["count_models"][0]

    @property
    def mdltyp(self) -> Optional[np.ndarray]:
        if not self.has_record(Record.MDLTYP):
            return
        return self.extract_record_fields(Record.MDLTYP)["description"]

    @property
    def author(self) -> Optional[np.ndarray]:
        if not self.has_record(Record.AUTHOR):
            return
        return self.extract_record_fields(Record.AUTHOR)["authors"]

    @property
    def revdat(self):
        if not self.has_record(Record.REVDAT):
            return
        fields = self.extract_record_fields(Record.REVDAT)
        unique_mod_nums, idx_unique_mod_nums = np.unique(
            fields["modification_num"], return_index=True
        )
        res_names = np.empty(unique_mod_nums.size, dtype=object)
        for i, unique_chain_id in enumerate(unique_mod_nums):
            mask_chain = fields["chain_id"] == unique_chain_id
            res_names_pre = fields["residue_name"][mask_chain][
                np.argsort(fields["serial"][mask_chain])
            ]
            res_names[i] = res_names_pre[res_names_pre != ""]
        df = pd.DataFrame(
            {
                "modification_type": fields["modification_type"][idx_unique_mod_nums],
                "residue_name": res_names
            }, index=unique_mod_nums
        )
        df.index.name = "chain_id"


    @property
    def sprsde(self) -> Optional[Tuple[datetime.date, np.ndarray]]:
        if not self.has_record(Record.SPRSDE):
            return
        fields = self.extract_record_fields(Record.SPRSDE)
        superseded_pdb_ids = fields["superseded_pdb_ids"][fields["superseded_pdb_ids"] != ""]
        return fields["superseded_date"][0], superseded_pdb_ids

    @property
    def jrnl(self):
        return

    @property
    def dbref(self) -> Optional[pd.DataFrame]:
        """
        Parse the DBREF and/or DBREF1/DBREF2 records of the PDB file.

        Returns
        -------
        pandas.DataFrame | None


        Notes
        -----
        Since both DBREF and DBREF1/DBREF2 records contain the same type of information,
        they are both parsed into the same dataframe.
        """
        if not (self.has_record(Record.DBREF) or self.has_record(Record.DBREF1)):
            return None
        # DBREF is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed. But information may be in both DBREF and
        # DBREF1/DBREF2 records. These are extracted separately, and merged:
        df_main = pd.DataFrame(columns=Record.DBREF.columns.keys())
        if self.has_record(Record.DBREF):
            # Extract DBREF records; no post-extraction processing needed.
            df_main = pd.concat(
                (df_main, self.extract_record_fields(Record.DBREF, to_dataframe=True))
            )
        if self.has_record(Record.DBREF1):
            # Extract DBREF1/DBREF2 records and merge rows to obtain the corresponding
            # DBREF-like dataframe.
            df = pd.concat(
                [
                    self.extract_record_fields(record, to_dataframe=True)
                    for record in (Record.DBREF1, Record.DBREF2)
                ],
                axis=1
            )
            df_main = pd.concat((df_main, df))
        return df_main

    @property
    def seqadv(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.SEQADV):
            return None
        # SEQADV is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed.
        return self.extract_record_fields(Record.SEQADV, to_dataframe=True)

    @property
    def seqres(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.SEQRES):
            return None
        # SEQRES is a multiple-times/multiple-lines record. Thus, after extraction, correlated
        # lines are merged to create a multiple-times/single-line dataframe format:
        fields = self.extract_record_fields(Record.SEQRES)
        # Continued lines all have the same chain ID
        unique_chain_ids, idx_unique_chain_ids = np.unique(fields["chain_id"], return_index=True)
        res_names = np.empty(unique_chain_ids.size, dtype=object)
        # Merge all lines for each chain:
        for i, unique_chain_id in enumerate(unique_chain_ids):
            mask_chain = fields["chain_id"] == unique_chain_id
            # Correlated lines are ordered according to the 'serial' field; this "should" always
            # start at 1 and increment by 1 each line, meaning the extracted lines should
            # already be in the correct order. Nevertheless, we re-order by the serial field,
            # just in case.
            res_names_pre = fields["residue_name"][mask_chain][
                np.argsort(fields["serial"][mask_chain])
            ]
            # Discard empty fields (these should be at the end of last line), and assign
            res_names[i] = res_names_pre[res_names_pre != ""]
        df = pd.DataFrame(
            {
                "residue_count": fields["residue_count"][idx_unique_chain_ids],
                "residue_name": res_names
            }, index=unique_chain_ids
        )
        df.index.name = "chain_id"
        return df

    def validate_seqres(self):
        # If no SEQRES record exists
        if not self.has_record(Record.SEQRES):
            # but ATOM records exist
            if self.has_record(Record.ATOM):
                self._raise_or_warn(
                    "SEQRES records not found, while ATOM records exist. "
                    "SEQRES is a mandatory record when ATOM records exist.",
                    raise_level=3
                )
        self._check_routine(Record.SEQRES)
        return

    @property
    def modres(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.MODRES):
            return
        # MODRES is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed.
        return self.extract_record_fields(Record.MODRES, to_dataframe=True)

    @property
    def het(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HET):
            return None
        # HET is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed.
        return self.extract_record_fields(Record.HET, to_dataframe=True)

    @property
    def hetnam(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HETNAM):
            return None
        # HETNAM is a multiple-times/multiple-lines record. Thus, after extraction, correlated
        # lines are merged to create a multiple-times/single-line dataframe format:
        fields = self.extract_record_fields(Record.HETNAM)
        _, idx_unique_ids = np.unique(fields["het_id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["het_id"][idx_unique_ids]
        names = list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["het_id"] == unique_id
            names_pre = fields["name"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            names.append(pdbfields.String.from_pdb(names_pre))
        df = pd.DataFrame(
            {
                "name": names
            }, index=unique_ids
        )
        df.index.name = "het_id"
        return df

    @property
    def hetsyn(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HETSYN):
            return None
        fields = self.extract_record_fields(Record.HETSYN)
        _, idx_unique_ids = np.unique(fields["het_id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["het_id"][idx_unique_ids]
        synonyms = list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["het_id"] == unique_id
            synonyms_pre = fields["synonyms"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            synonyms.append(pdbfields.SList.from_pdb(pdbfields.String.from_pdb(synonyms_pre)))
        df = pd.DataFrame(
            {
                "synonyms": synonyms
            }, index=unique_ids
        )
        df.index.name = "het_id"
        return df

    @property
    def formul(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.FORMUL):
            return None
        fields = self.extract_record_fields(Record.FORMUL)
        _, idx_unique_ids = np.unique(fields["het_id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["het_id"][idx_unique_ids]
        is_water, formula, num_in_chain, num_out_chain = list(), list(), list(), list()
        for i, unique_id in enumerate(unique_ids):
            mask_chain = fields["het_id"] == unique_id
            is_water.append(np.any(fields["is_water"][mask_chain]))
            formula_pre = fields["formula"][mask_chain][
                np.argsort(fields["continuation"][mask_chain])
            ]
            # All occurrences of the HET group within a chain are grouped together with a
            # multiplier. The remaining occurrences are also grouped with a multiplier, e.g.:
            # 4(2(C21 H28 N7 O17 P3)), 2(C19 H19 N7 O6) , or C10 H22 O6.
            formula_single_line = pdbfields.String.from_pdb(formula_pre)
            formula_split = re.split(r'[()]+', formula_single_line)
            num_split = len(formula_split)
            if num_split == 1:
                formula.append(formula_split[0])
                num_in_chain.append(0)
                num_out_chain.append(0)
            elif num_split == 3:
                formula.append(formula_split[-2])
                num_in_chain.append(formula_split[0])
                num_out_chain.append(0)
            elif num_split == 4:
                formula.append(formula_split[-2])
                num_in_chain.append(formula_split[1])
                num_out_chain.append(formula_split[0])
            else:
                raise PDBParsingError(
                    f"Could not parse FORMUL record's field {formula_single_line}."
                )
        df = pd.DataFrame(
            {
                "component_num": fields["component_num"][idx_unique_ids],
                "is_water": is_water,
                "formula": formula,
                "count_in_chain": num_in_chain,
                "count_outer_chain": num_out_chain,
            }, index=unique_ids
        )
        df.index.name = "het_id"
        return df

    @property
    def helix(self) -> Optional[pd.DataFrame]:
        if not self.has_record(Record.HELIX):
            return
        # HELIX is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed.
        return self.extract_record_fields(Record.HELIX, to_dataframe=True)

    @property
    def sheet(self):
        if not self.has_record(Record.SHEET):
            return
        # SHEET is a multiple-times/one-line record; i.e. each line is independent of the others.
        # Thus, no post-extraction processing is needed.
        return self.extract_record_fields(Record.SHEET, to_dataframe=True)

    @property
    def ssbond(self):
        if not self.has_record(Record.SSBOND):
            return
        return self.extract_record_fields(Record.SSBOND, to_dataframe=True)

    @property
    def link(self):
        return

    @property
    def cispep(self):
        return

    @property
    def site(self):
        if not self.has_record(Record.SITE):
            return
        fields = self.extract_record_fields(Record.SITE)
        _, idx_unique_ids = np.unique(fields["site_id"], return_index=True)
        idx_unique_ids = np.sort(idx_unique_ids)
        unique_ids = fields["site_id"][idx_unique_ids]
        site_ids, res_counts, res_names, chain_ids, res_nums, res_icodes = [], [], [], [], [], []
        for i, unique_id in enumerate(unique_ids):
            mask_site = fields["site_id"] == unique_id
            res_names_site = fields["residue_name"][mask_site]
            mask_empty = res_names_site != ""
            res_names_site = res_names_site[mask_empty]
            num_residues = res_names_site.size
            site_ids.extend([unique_id] * num_residues)
            res_counts.extend([fields["residue_count"][mask_site][0]] * num_residues)
            res_names.append(res_names_site)
            chain_ids.append(fields["chain_id"][mask_site][mask_empty])
            res_nums.append(fields["residue_num"][mask_site][mask_empty])
            res_icodes.append(fields["residue_icode"][mask_site][mask_empty])
        df = pd.DataFrame(
            {
                "residue_count": res_counts,
                "residue_name": np.concatenate(res_names),
                "chain_id": np.concatenate(chain_ids),
                "residue_num": np.concatenate(res_nums),
                "residue_icode": np.concatenate(res_icodes),
            }, index=site_ids
        )
        df.index.name = "site_id"
        return df

    @property
    def cryst1(self):
        if not self.has_record(Record.CRYST1):
            return
        fields = self.extract_record_fields(Record.CRYST1)
        unit_cell_lengths = np.array([fields["a"][0], fields["b"][0], fields["c"][0]])
        unit_cell_angles = np.array([fields["alpha"][0], fields["beta"][0], fields["gamma"][0]])
        return records.Cryst1(
            unit_cell_lengths=unit_cell_lengths,
            unit_cell_angles=unit_cell_angles,
            space_group=fields["space_group"][0],
            z=fields["z"][0]
        )

    @property
    def origx(self):
        if not self.has_record(Record.ORIGX1):
            return
        return self._parse_origx_scale("origx")

    @property
    def scale(self):
        if not self.has_record(Record.SCALE1):
            return
        return self._parse_origx_scale("scale")

    def _parse_origx_scale(self, record: Literal["origx", "scale"]):
        records = (
            Record.ORIGX1, Record.ORIGX2, Record.ORIGX3
        ) if record == "origx" else (Record.SCALE1, Record.SCALE2, Record.SCALE3)
        l1, l2, l3 = (
            self.record_lines(rec, as_char=False)[0]
            for rec in records
        )
        transform_matrix = np.array([[x[10:20], x[20:30], x[30:40]] for x in (l1, l2, l3)]).astype(np.float_)
        translation_vector = np.array([x[45:55] for x in (l1, l2, l3)]).astype(np.float_)
        return transform_matrix, translation_vector

    @property
    def mtrix(self):
        return

    @property
    def atom(self):
        return

    @property
    def anisou(self) -> Optional[pd.DataFrame]:
        return self._parse_record_in_model(Record.ANISOU)

    @property
    def ter(self) -> Optional[pd.DataFrame]:
        return self._parse_record_in_model(Record.TER)

    @property
    def hetatm(self):
        return

    @property
    def conect(self):
        return

    @property
    def master(self):
        if not self.has_record(Record.MASTER):
            return
        line = self.record_lines(Record.MASTER, as_char=False)[0]
        for start, stop in ((10, 15), (20, 25), (25, 30), (30, 35)):
            int(line[start:stop])

        return

    def _parse_record_in_model(self, record: Record) -> Optional[pd.DataFrame]:
        if not self.has_record(record):
            return
        fields = self.extract_record_fields(record)
        if not self.has_record(Record.MODEL):
            model_nr = np.full(shape=self.count_record(record), fill_value=1, dtype=np.ubyte)
        else:
            count_record_per_model = [
                np.count_nonzero(
                    np.logical_and(
                        self.indices_record_lines(record) > model_start_idx,
                        self.indices_record_lines(record) < model_end_idx,
                    )
                ) for model_start_idx, model_end_idx in zip(
                    self.indices_record_lines(Record.MODEL), self.indices_record_lines(Record.ENDMDL)
                )
            ]
            model_nr = np.repeat(
                np.arange(1, self.count_record(Record.MODEL) + 1), repeats=count_record_per_model
            )
        return pd.DataFrame({"model_nr": model_nr, **fields})

    def validate_master(self) -> sections.Bookkeeping:
        """
        Parse the MASTER record in a PDB file.

        Returns
        -------
        BookkeepingSection

        Raises
        ------
        PDBParsingError

        Notes
        -----
        MASTER record is a mandatory one-time/single-line record that must appear only once, as the
        second-last record in every PDB file. It is a control record, counting the occurrence
        (checksum) of selected record types in the PDB file.
        """
        # Count relevant records to construct the MASTER record data
        bookkeeping_section = pdbstruct.BookkeepingSection(
            remark=np.count_nonzero(self._lines_rec_names == "REMARK"),
            het=np.count_nonzero(self._lines_rec_names == "REMARK"),
            helix=np.count_nonzero(self._lines_rec_names == "HELIX"),
            sheet=np.count_nonzero(self._lines_rec_names == "SHEET"),
            site=np.count_nonzero(self._lines_rec_names == "SITE"),
            xform=np.count_nonzero(
                np.isin(
                    self._lines_rec_names,
                    [
                        "ORIGX1",
                        "ORIGX2",
                        "ORIGX3",
                        "SCALE1",
                        "SCALE2",
                        "SCALE3",
                        "MTRIX1",
                        "MTRIX2",
                        "MTRIX3",
                    ]
                )
            ),
            coord=np.count_nonzero(np.isin(self._lines_rec_names, ["ATOM", "HETATM"])),
            ter=np.count_nonzero(self._lines_rec_names == "TER"),
            connect=np.count_nonzero(self._lines_rec_names == "CONECT"),
            seqres=np.count_nonzero(self._lines_rec_names == "SEQRES"),
        )
        # Look for existing MASTER records in file
        idx_master_record = np.argwhere(self._lines_rec_names == "MASTER").rehape(-1)
        # Case 1: no existing MASTER record in file
        if idx_master_record.size == 0:
            self._raise_or_warn(
                "MASTER record not found. It is a mandatory record and must appear once in "
                "the file, as the second-last record (before the END record).",
                raise_level=3,
            )
            return bookkeeping_section

        if idx_master_record.size > 1:
            self._raise_or_warn(
                f"Multiple MASTER records found ({idx_master_record.size} at lines "
                f"{idx_master_record + 1}). It is a one-time/single-line record and must appear "
                "only once in the file.",
                raise_level=2,
            )
        # Parse existing MASTER record
        master = parsing.extract_columns(self._lines_chars[idx_master_record[-1]],
                                         Record.MASTER.columns)
        master_record = pdbstruct.BookkeepingSection(*master)

        if master_record != bookkeeping_section:
            self._raise_or_warn(
                "MASTER record values do not match with number of records in file. "
                f"MASTER record values are {master_record}, "
                f"but actual record counts are {bookkeeping_section}.",
                raise_level=2,
            )

        if idx_master_record != self._lines_rec_names.size - 2:
            self._raise_or_warn(
                f"MASTER record found in wrong position (line {idx_master_record[0] + 1}). "
                f"It must be the second-last record in the file, "
                f"i.e. at line {self._lines_rec_names.size - 2}.",
                raise_level=3,
            )
        empty_cols = np.concatenate((np.arange(7, 11), np.arange(71, 81))) - 1
        if np.any(self._lines_chars[idx_master_record[0], empty_cols] != " "):
            self._raise_or_warn(
                "MASTER record line contains characters in columns that must be empty (i.e. "
                f"columns 7–10 and 71–80). "
                f"Found: {self._lines_chars[idx_master_record[0], empty_cols]}",
                raise_level=3,
            )
        return bookkeeping_section

    def validate_section_coordinate(self):

        if ind_ter.size == 0:
            if ind_atom.size != 0:
                self._raise_or_warn(
                    "ATOM records found without TER record.",
                    raise_level=1,
                )

        indices_record_model = self._idx_record[Record.MODEL]  # line indices of MODEL records
        count_model = indices_record_model.size  # count of MODEL records
        if count_model > 0:
            indices_record_endmdl = self._idx_record[Record.ENDMDL]
            count_endmdl = indices_record_endmdl.size
            if count_model != count_endmdl:
                self._raise_or_warn(
                    "Number of MODEL and ENDMDL records do not match: "
                    f"{count_model} vs {count_endmdl}",
                    raise_level=1
                )
            model_lengths = indices_record_endmdl - indices_record_model
            if model_lengths[0] < 1 or np.any(model_lengths != model_lengths[0]):
                self._raise_or_warn(
                    "Models have different lengths or are empty",
                    raise_level=1
                )

        else:
            # If there is no MODEL record, then there is only one model
            # Extract the records of the single model
            filter_model = np.isin(self._lines_rec_names, ["ATOM", "ANISOU", "TER", "HETATM"])

            model = read_model(rec_names[filter_model], rec_chars[filter_model])
            return CoordinateSection(models=(model,))

        models = []
        for start_idx, end_idx in zip(models_start_ind, models_end_ind):
            models.append(
                read_model(
                    record_types=rec_names[start_idx + 1:end_idx],
                    pdb_lines=rec_chars[start_idx + 1:end_idx]
                )
            )
        return CoordinateSection(models=tuple(models))


    def validate_dbref(self):
        if not (
                self._has_record[Record.DBREF]
                or self._has_record[Record.DBREF1]
                or self._has_record[Record.DBREF2]
        ):
            if self._has_record[Record.ATOM] or self._has_record[Record.TER]:
                self._raise_or_warn(
                    "DBREF and DBREF1/DBREF2 records not found, while ATOM records exist. "
                    "At least one must exist when ATOM records exist.",
                    raise_level=3
                )
        if self._has_record[Record.DBREF1]:
            if self._idx_record[Record.DBREF1].size != self._idx_record[Record.DBREF2].size:
                self._raise_or_warn(
                    "Number of DBREF1 and DBREF2 records do not match. DBREF1/DBREF2 records are "
                    "two-line records that always have to be present together.",
                    raise_level=2
                )
            self.validate_empty_columns(Record.DBREF1)
            self.validate_empty_columns(Record.DBREF2)
        self._check_routine(Record.DBREF)
        return

    def validate_end(self) -> NoReturn:
        """
        Validate the END record in the PDB file.

        Raises
        ------
        PDBParsingError

        Notes
        -----
        END record is a mandatory one-time/single-line record that must appear only once,
        as the final record at the end of every PDB file.
        The first 3 characters of the line must be 'END', followed by 77 spaces to make
        an 80-character line. Since it contains no information, it can be completely ignored when
        PDB file-format validation is not a concern.
        """
        if not self.has_record(Record.END):
            self._raise_or_warn(
                "END record not found. It is a mandatory record and must appear once at "
                "the end of the file.",
                raise_level=3,
            )
        idx_record_end = self.indices_record_lines(Record.END)
        if self.count_record(Record.END) > 1:
            self._raise_or_warn(
                f"Multiple END records found ({idx_record_end.size} at lines "
                f"{idx_record_end + 1}). It is a one-time/single-line record and must appear "
                "only once in the file.",
                raise_level=3,
            )
        if idx_record_end != self.count_lines - 1:
            self._raise_or_warn(
                f"END record found at the middle of the file (line {idx_record_end[0] + 1}). "
                "It must be the final record in the file.",
                raise_level=3,
            )
        lines_record_end = self.record_lines(Record.END, drop_name=True)
        mask_non_empty = lines_record_end != " "
        if np.any(mask_non_empty):
            self._raise_or_warn(
                "END record line contains extra characters. It must contain only the record name "
                "'END' at beginning of the line. "
                f"Found characters {lines_record_end[mask_non_empty]} at positions "
                f"{np.argwhere(mask_non_empty)}.",
                raise_level=3,
            )
        return

    def _validate_existence_mandatory_record(self, record: Record):
        if not self.has_record(record):
            self._raise_or_warn(
                f"{record.name} record not found. It is a mandatory record.",
                raise_level=3
            )
            return False
        return True

    def _validate_one_time_single_line_record(self, record: Record) -> bool:
        if self.count_record(record) > 1:
            self._raise_or_warn(
                f"Multiple {record.name} records found at lines "
                f"{', '.join(self.indices_record_lines(record))}. "
                f"{record.name} is a one-time/single-line record and must only appear once in the file.",
                raise_level=1
            )
            return False
        return True

    def _validate_pdb_ids(self, pdb_ids: np.ndarray):
        validated = True
        for pdb_id in pdb_ids:
            if not pdb_id[0].isdigit():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. The first character in the ID must be a digit.",
                    raise_level=3
                )
                validated = False
            if not pdb_id[1:].isalpha():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. The last three characters in the ID must be alphabetical.",
                    raise_level=3
                )
                validated = False
            if not pdb_id.isupper():
                self._raise_or_warn(
                    f"PDB ID {pdb_id} not correct. All characters in the ID must be upper-case.",
                    raise_level=3
                )
                validated = False
        return validated





    def validate_empty_columns(self, record: Record):
        if self.has_record(record):
            non_spaces = self.record_lines(record)[:, record.empty_columns] != " "
            if np.any(non_spaces):
                positions = ", ".join([f"({x}, {y})" for x, y in np.argwhere(non_spaces)])
                self._raise_or_warn(
                    f"{record.name} records contain non-space characters in columns that must be empty."
                    f"Positions (line number, character number): {positions}",
                    raise_level=1
                )
                return False
        return True

    def validate_line_width(self):
        # Check if all lines are 80 characters:
        line_widths = np.char.str_len(self._lines)
        lines_not80 = line_widths != 80
        if np.any(lines_not80):
            faulty_lines = [
                f"{line_num} ({line_width})" for line_num, line_width in zip(
                    np.argwhere(lines_not80).reshape(-1) + 1,
                    line_widths[lines_not80]
                )
            ]
            self._raise_or_warn(
                "Following lines are not 80 characters long (lengths given in brackets): "
                f"{', '.join(faulty_lines)}",
                raise_level=3
            )
        return

    def validate_record_names(self):
        if not np.all(np.isin(self._lines_rec_names, Record.names)):
            self._raise_or_warn(
                f"Unrecognized record names found.",
                raise_level=3
            )

    def mask_record(self, record: Record) -> np.ndarray:
        """
        A boolean array indicating whether each line in the PDB file is of a given record type.

        Parameters
        ----------
        record

        Returns
        -------
        ndarray
        """
        mask = self._mask_record.get(record)
        if mask is not None:
            return mask
        mask = self._lines_rec_names == record.name
        self._mask_record[record] = mask
        return mask

    def indices_record_lines(self, record: Record) -> np.ndarray:
        idx = self._idx_record.get(record)
        if idx is not None:
            return idx
        idx = np.argwhere(self.mask_record(record)).reshape(-1)
        self._idx_record[record] = idx
        return idx

    def count_record(self, record: Record) -> int:
        count = self._count_record.get(record)
        if count is not None:
            return count
        count = self.indices_record_lines(record).size
        self._count_record[record] = count
        return count

    def has_record(self, record: Record) -> bool:
        has = self._has_record.get(record)
        if has is not None:
            return has
        has = self.count_record(record) != 0
        self._has_record[record] = has
        return has

    def record_lines(self, record: Record = None, as_char: bool = True, drop_name: bool = False):
        idx = self.indices_record_lines(record) if record is not None else slice(None)
        if as_char:
            return self._lines_chars[idx, (6 if drop_name else 0):]
        return self._lines_rec_content[idx] if drop_name else self._lines[idx]

    def extract_record_fields(
            self,
            record: Record,
            to_dataframe: bool = False
    ) -> Union[Dict[str, Any], pd.DataFrame]:
        lines = self.record_lines(record=record)
        return parsing.extract_columns(
            char_array=lines, columns=record.columns, to_dataframe=to_dataframe
        )

    @property
    def count_lines(self) -> int:
        """Total number of lines in the PDB file."""
        return self._lines.size

    @property
    def strictness(self) -> Literal[0, 1, 2, 3]:
        """
        Level of strictness for raising exceptions and warnings when encountering mistakes
        in the PDB file:
            0: raise only fatal errors and don't show any warnings.
            1: raise only fatal errors. All other errors are reported as warnings.
            2: raise fatal errors and mistakes resulting in ambiguous data. Inconsequential mistakes are
            reported as warnings.
            3: completely validate the PDB file and raise all errors.

        Returns
        -------
        literal[0, 1, 2, 3]
        """
        return self._strictness

    @strictness.setter
    def strictness(self, value):
        if value not in (0, 1, 2, 3):
            raise ValueError(
                "Parameter `strictness` expects a value of either 0, 1, 2, or 3. "
                f"Input argument was: {value}."
            )
        self._strictness = value
        return

    def _check_routine(self, record: Record):
        if self.has_record(record):
            self._check_records_are_consecutive(record)
            self.validate_empty_columns(record)



    def _check_records_are_consecutive(self, record: Record):
        if self._has_record[record] and np.any(np.diff(self._idx_record[record]) != 1):
            self._raise_or_warn(
                f"{record.name} records are not in consecutive lines.",
                raise_level=3
            )
        return

    def _raise_or_warn(
            self,
            msg: str,
            raise_level: Literal[1, 2, 3]
    ) -> NoReturn:
        # if self.strictness >= raise_level:
        #     raise PDBParsingError(msg)
        # else:
        warnings.warn(msg)
        return


def from_pdb_id(
        pdb_id: str,
        biological_assembly_id: Optional[int] = None,
        validate: bool = False
) -> PDBParser:
    """
    Create a PDB parser from a given PDB ID.


    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    biological_assembly_id : int | str, optional, default: None
        Biological assembly ID. If not provided (i.e. when `None`), the asymmetric unit will be
        downloaded. Notice that many records are only available in the PDB file of the asymmetric unit.
    validate

    Returns
    -------
    opencadd.io.pdb.Parser
    """
    pdb_file = oc.api.web.rcsb.file_pdb_entry(
        pdb_id=pdb_id, file_format="pdb", biological_assembly_id=biological_assembly_id
    )
    return PDBParser(file=pdb_file, strictness=validate)





# pdbfile = rcsb.file_pdb_entry("3w32", "pdb")
# parser = PDBParser(pdbfile)#("/Users/home/Downloads/3w32as.pdb")
# x= parser.compnd
# print(x, type(x))