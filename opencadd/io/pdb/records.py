

from typing import Tuple, Sequence
import datetime

import numpy as np
import pandas as pd

from opencadd.io.parsing import extract_column
from opencadd.io.pdb import fields


class Header:
    """
    HEADER record contains the PDB ID, classification, and deposition date of the entry.

    Notes
    -----
    HEADER is:
        1. mandatory; it must be present in all PDB files.
        2. one-time/single-line: it only appears once in the PDB file.
        3. fixed: must always be the first record in the file.

    References
    ----------
    https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#HEADER
    """
    def __init__(
            self,
            classification: Tuple[Tuple[str]],
            deposition_date: datetime.date,
            pdb_id: str,
    ):
        self._classification: Tuple[Tuple[str]] = classification
        self._deposition_date: datetime.date = deposition_date
        self._pdb_id: str = pdb_id
        return

    def __repr__(self):
        return f"HEADER({self.classification}, {self.deposition_date}, {self.pdb_id})"

    def __str__(self):
        return f"PDB ID {self.pdb_id} (deposited on {self.deposition_date}): {self.classification}"

    @property
    def pdb_id(self) -> str:
        """
        PDB ID of the entry; this identifier is unique within the PDB.
        """
        return self._pdb_id

    @property
    def classification(self) -> Tuple[Tuple[str]]:
        """
        Classification of each molecule within the entry.

        Returns
        -------
        Tuple[Tuple[str]]
            Each sub-tuple corresponds to one molecule in the entry, with each element
            describing one classification/function.

        Notes
        -----
        * Classification may be based on function, metabolic role, molecule type, cellular location
        etc., but it exactly matches one of a collection of strings, available at:
        https://www.wwpdb.org/docs.html
        *  Due to the limited length of the classification field, strings must sometimes be
        abbreviated. In these cases, the full terms are given in KEYWDS records in no strict order.
        """
        return self._classification

    @property
    def deposition_date(self) -> datetime.date:
        """
        Date of deposition of the entry at the Protein Data Bank, i.e. the 'depDate' field.
        """
        return self._deposition_date

    @property
    def pdb_format(self) -> str:
        classification = "/".join(", ".join(entry) for entry in self.classification)
        dep_date = fields.Date.to_pdb(self.deposition_date)
        return f"HEADER{'':4}{classification:<40}{dep_date}{'':3}{self.pdb_id}{'':14}"


class Obsolete:
    def __init__(
            self,
            pdb_id: str,
            replacement_date: datetime.date,
            replaced_pdb_ids: Sequence[str],
    ):
        self.pdb_id: str = pdb_id
        self.replacement_date: datetime.date = replacement_date
        self.replaced_pdb_ids: np.ndarray = replaced_pdb_ids
        return




class Split:
    pass


class Caveat:
    pass



class Source:
    pass


class Keywords:
    pass


class ExperimentalData:
    pass


class NumModels:
    pass


class ModelType:
    pass


class Author:
    pass


class RevisionData:
    pass


class SPRSDE:
    pass


class Journal:
    pass


class Remark2:
    pass


class Remark3:
    pass
