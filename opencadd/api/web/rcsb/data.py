"""
Functions for retrieving data from the PDB archive.

References
----------
https://data.rcsb.org/
https://data.rcsb.org/redoc/
https://data.rcsb.org/data-attributes.html
"""

from typing import Literal

import numpy as np
from opencadd.api.web.web import response_http_request

_URL_ROOT = "https://data.rcsb.org/rest/v1/"
_URL_HOLDINGS = f"{_URL_ROOT}holdings/"
_URL_HOLDINGS_CURRENT_PDB_IDS = f"{_URL_HOLDINGS}current"


class APIConst:
    class Data:
        # Prefix for all data APIs
        URL_PREFIX = "https://data.rcsb.org/rest/v1/"

        class Entry:
            """
            Provides access to information about PDB structures at the
            top level of PDB structure hierarchical data organization.
            """

            # Get PDB structure data by PDB-ID
            URL_STRUCTURE = "core/entry/"
            # i.e. core/entry/{entry_id}

            # Get PubMed annotations (data integrated from PubMed)
            # for a given entry's primary citation.
            URL_PUBMED = "core/pubmed/"
            # i.e. core/pubmed/{entry_id}

        class Entity:
            """
            Provides access to information about PDB structures
            at the level of unique macromolecular entities.
            """

            # Get polymer entity data by PDB ID and ENTITY ID.
            URL_POLYMER = "core/polymer_entity/"
            # i.e. core/polymer_entity/{entry_id}/{entity_id}

            # Get branched entity description by PDB ID and ENTITY ID.
            URL_BRANCHED = "core/branched_entity/"
            # i.e. core/branched_entity/{entry_id}/{entity_id}

            # Get non-polymer entity data by PDB ID and ENTITY ID.
            URL_NON_POLYMER = "core/nonpolymer_entity/"
            # i.e. core/nonpolymer_entity/{entry_id}/{entity_id}

            # Get UniProt annotations for a given macromolecular
            # entity (identified by PDB ID and ENTITY ID).
            URL_UNIPROT = "core/uniprot/"
            # i.e. core/uniprot/{entry_id}/{entity_id}

        class EntityInstance:
            """
            Provides access to information about PDB structures
            at the level of unique molecular entities (e.g. chains).
            """

            # Get polymer entity instance (a.k.a chain) data
            # by PDB ID and ASYM ID (label_asym_id).
            URL_POLYMER = "core/polymer_entity_instance/"
            # i.e. core/polymer_entity_instance/{entry_id}/{asym_id}

            # Get branched entity instance description by PDB ID and ASYM ID.
            URL_BRANCHED = "core/branched_entity_instance/"
            # i.e. core/branched_entity_instance/{entry_id}/{asym_id}

            # Get non-polymer entity instance description
            # by PDB ID and ASYM ID (label_asym_id).
            URL_NON_POLYMER = "core/nonpolymer_entity_instance/"
            # i.e. core/nonpolymer_entity_instance/{entry_id}/{asym_id}

        class Assembly:
            """
            Provides access to information about PDB structures
            at the quaternary structure level.
            """

            # Get structural assembly description by PDB ID and ASSEMBLY ID.
            URL_ASSEMBLY = "core/assembly/"
            # i.e. core/assembly/{entry_id}/{assembly_id}

        class Chemical:
            """
            Provides access to information about chemical components
            from which the relevant chemical structures can be constructed.
            """

            # Get chemical component (ligands, small molecules and monomers)
            # by CHEM COMP ID (CCD ID) that uniquely identifies the chemical
            # component. For protein polymer entities, this is the three-letter
            # code for the amino acid. For nucleic acid polymer entities,
            # this is the one-letter code for the base.
            URL_COMPONENT = "core/chemcomp/"
            # i.e. core/chemcomp/{comp_id}

            # Get DrugBank annotations (integrated from DrugBank resource)
            # for a given chemical component (identified by CHEM COMP ID).
            URL_DRUG_BANK = "core/drugbank/"
            # i.e. core/drugbank/{comp_id}

    class File:
        PDB = "https://files.rcsb.org/download/"
        SMALL_MOLECULE = "https://files.rcsb.org/ligands/download/"


def get_data_entry(data_type, pdb_code):
    """
    Provides information about PDB structures at the top
    level of PDB structure hierarchical data organization.
    """
    return _get_data("Entry", data_type, pdb_code)


def get_data_entity(data_type, pdb_code, entity_id):
    return _get_data("Entity", data_type, pdb_code, entity_id)


def get_data_entity_instance(data_type, pdb_code, instance_id):
    return _get_data("EntityInstance", data_type, pdb_code, instnace_id)


def get_data_assembly(pdb_code, assembly_id):
    return _get_data("Assembly", "assembly", pdb_code, assembly_id)


def get_data_chemical(data_type, compound_id):
    return _get_data("Chemical", data_type, compound_id)


def _get_data(data_class, data_type, *args):
    full_url = _build_data_url(data_class, data_type, *args)
    data = _send_request(full_url)
    return data


def _build_data_url(data_class, data_type, *args):
    url_prefix = APIConst.Data.URL_PREFIX
    url = getattr(getattr(APIConst.Data, data_class), f"URL_{data_type.upper()}")
    url_suffix = "".join([f"{value}/" for value in args])
    full_url = f"{url_prefix}{url}{url_suffix}"
    return full_url


def pdb_id_holdings(data: Literal["current", "unreleased", "removed"] = "current") -> np.ndarray:
    """
    Get all PDB ID holdings data.

    Parameters
    ----------
    data : Literal["current", "unreleased", "removed"], optional, default: "current"
        Data to retrieve; either all current, all unreleased, or all removed PDB IDs.

    Returns
    -------
    numpy.ndarray[ndim=1, dtype="<U4"]
        A one-dimensional array of PDB IDs.
    """
    if data not in ("current", "unreleased", "removed"):
        raise ValueError(f"{data} is not a valid argument for `data`.")
    return np.array(
        response_http_request(url=f"{_URL_HOLDINGS}{data}/entry_ids", response_type="json")
    )
