from typing import NamedTuple, List, Tuple, Dict, Sequence
import datetime
import warnings

import numpy as np
import pandas as pd


class HeaderRecord(NamedTuple):
    pdb_id: str = None
    classification: Tuple[Tuple[str]] = None
    deposition_date: datetime.date = None
    invalid_lines: Tuple = None






class SeqResRecord:
    COL_RANGES = np.concatenate(
        (
            np.array([[8, 10], [12, 12], [14, 17]]),
            np.stack((np.arange(20, 69, 4), np.arange(22, 71, 4)), axis=-1)
        )
    )

    def __init__(self, sequences: Dict[str, Sequence[str]]):
        self._data = sequences
        return

    @property
    def chain_ids(self) -> List:
        pass

    @classmethod
    def from_character_table(cls, record_types: np.ndarray, char_tab: np.ndarray):
        is_record_line: np.ndarray = record_types == "SEQRES"
        record_lines = char_tab[is_record_line]

        # Get 'serNum' column
        serial_nums = record_lines[:, 7:10].view(dtype=(str, 3)).reshape(-1).astype(int)
        # Get 'chainID' column
        chain_ids = record_lines[:, 11:12].reshape(-1)
        # Get 'numRes' column
        num_res = record_lines[:, 13:17].view(dtype=(str, 4)).reshape(-1).astype(int)



class TitleSection(NamedTuple):
    header: HeaderRecord


class RemarkSection(NamedTuple):
    pass


class PrimaryStructureSection(NamedTuple):
    pass


class HeterogenSection(NamedTuple):
    pass


class SecondaryStructureSection(NamedTuple):
    pass


class ConnectivityAnnotationSection(NamedTuple):
    pass


class SiteSection(NamedTuple):
    pass


class CrystallographicSection(NamedTuple):
    pass


class CoordinateTransformationSection(NamedTuple):
    pass


class CoordinateSection(NamedTuple):

    atom: pd.DataFrame = None
    anisou: pd.DataFrame = None
    hetatm: pd.DataFrame = None
    ter: pd.DataFrame = None


class ConnectivitySection(NamedTuple):
    pass


class BookkeepingSection(NamedTuple):
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


class PDBFile(NamedTuple):
    title: TitleSection
    remark: RemarkSection
    primary_structure: PrimaryStructureSection
    heterogen: HeterogenSection
    secondary_structure: SecondaryStructureSection
    connectivity: ConnectivitySection
    site: SiteSection
    crystallographic: CrystallographicSection
    coordinate: CoordinateSection
    bookkeeping: BookkeepingSection

