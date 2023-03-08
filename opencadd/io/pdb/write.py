from typing import Union, Sequence, Optional
from pathlib import Path
import pandas as pd
from opencadd import _typing, _sysio
import numpy as np
import opencadd as oc
from . import struct, _writer


def from_chemsys(
        system,
        output_filename: Optional[str] = None,
        output_path: _typing.PathLike = None,
        models: Optional[Union[int, Sequence[int]]] = None,
        separate_models: bool = False,
):
    pdb_struct = struct.PDBStructure(
        atom=system.composition
    )
    return from_structure(
        pdb_structure=pdb_struct,
        output_filename=output_filename,
        output_path=output_path,
        models=models,
        separate_models=separate_models
    )


def from_structure(
        pdb_structure: struct.PDBStructure,
        output_filename: Optional[str] = None,
        output_path: _typing.PathLike = None,
        models: Optional[Union[int, Sequence[int]]] = None,
        separate_models: bool = False,
):
    writer = _writer.PDBWriter(pdb_structure=pdb_structure)
    return writer.write(
        models=models,
        separate_models=separate_models,
    )


# output_filepaths = []
#             for model_num, pdb_str in zip(selected_models, pdb_strings):
#                 output_filepaths.append(
#                     _sysio.save_to_file(
#                         content=pdb_str,
#                         filename=f"{output_filename}_{model_num}",
#                         extension="pdb",
#                         path=output_path
#                     )
#                 )
#             return output_filepaths