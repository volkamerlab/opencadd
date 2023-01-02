


def write_header_record(
        header: HeaderRecord,
        valid_only: bool = True,
        default_pdb_id: str = "XXXX",
        default_classification: str = "N/A",
        default_deposition_date: str = "??-???-??",
) -> str:
    pdb_id = header.pdb_id if header.pdb_id is not None else default_pdb_id
    classification = "/".join(
        ", ".join(entity_functions) for entity_functions in header.classification
    ) if header.classification is not None else default_classification
    deposition_date = header.deposition_date.strftime("%d-%b-%y").upper(
    ) if header.deposition_date is not None else default_deposition_date
    valid_header = f"HEADER{'':4}{classification:<40}{deposition_date}{'':3}{pdb_id}{'':14}"
    return valid_header if valid_only else valid_header + "\n" + "\n".join(header.invalid_lines)


def record_master() -> str:
    pass


def record_end() -> str:
    """
    END record of a PDB file.
    The END record marks the end of the PDB file, and must appear as the final record in every
    file.
    """
    return f"{'END':<80}"
