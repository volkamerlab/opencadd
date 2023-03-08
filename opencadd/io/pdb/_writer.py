from typing import Optional, Union, Sequence

import pandas as pd
import numpy as np

import opencadd as oc
from opencadd import _typing, _sysio
from . import struct


__author__ = "Armin Ariamajd"


class PDBWriter:
    def __init__(
            self,
            pdb_structure: struct.PDBStructure
    ):
        self._struct = pdb_structure
        self._pre_coord_records_str = None
        self._models_str = None
        self._post_coord_records_str = None
        return

    def write(self, models, separate_models):

        if self._pre_coord_records_str is None:
            self._pre_coord_records_str = self.pre_coord_str()
        if self._post_coord_records_str is None:
            self._post_coord_records_str = self.post_coord_str()

        if models is None:
            selected_models = self._struct.atom.model_num.unique()
        elif np.issubdtype(type(models), np.integer):
            selected_models = [models]
        else:
            selected_models = models

        models_str = self.section_coordinate(models=selected_models)
        if separate_models:
            pdb_strings = [
                f"{self._pre_coord_records_str}{model_str}{self._post_coord_records_str}"
                for model_str in models_str
            ]
            return pdb_strings

        coordinate_section = "\n".join(models_str)
        pdb_string = f"{self._pre_coord_records_str}{coordinate_section}{self._post_coord_records_str}"
        return pdb_string

    def pre_coord_str(self):
        return ""

    def post_coord_str(self):
        return f"\n{'END':<80}"

    def section_coordinate(
            self,
            models: Sequence[int],
    ):
        models_str = []
        atom_df = self._struct.atom
        for model_num in models:
            pdb_lines = [f"MODEL{'':5}{model_num:>4}{'':66}"]
            model = atom_df[atom_df.model_num == model_num]
            for chain_id in model.chain_id.unique():
                chain = model[model.chain_id == chain_id]
                poly = chain[chain.res_poly]
                for i in poly.index:
                    row = poly.loc[i]
                    pdb_lines.append(
                        self.line_coordinates(
                            is_std=row.res_std,
                            serial=row.serial,
                            atom_name=row.atom_name,
                            alt_loc=row.alt_loc,
                            res_name=row.res_name,
                            chain_id=chain_id,
                            res_num=row.res_num,
                            res_icode=row.res_icode,
                            x=row.x,
                            y=row.y,
                            z=row.z,
                            occupancy=row.occupancy,
                            temp_factor=row.temp_factor,
                            element=row.element,
                            charge=row.charge
                        )
                    )
                pdb_lines.append(
                    f"TER   {row.serial + 1:>5}{'':6}{row.res_name:>3} {chain_id:1}"
                    f"{row.res_num:>4}{row.res_icode:1}{'':53}"
                )
            hets = model[~model.res_poly]
            for i in hets.index:
                row = hets.loc[i]
                pdb_lines.append(
                    self.line_coordinates(
                        is_std=False,
                        serial=row.serial,
                        atom_name=row.atom_name,
                        alt_loc=row.alt_loc,
                        res_name=row.res_name,
                        chain_id=chain_id,
                        res_num=row.res_num,
                        res_icode=row.res_icode,
                        x=row.x,
                        y=row.y,
                        z=row.z,
                        occupancy=row.occupancy,
                        temp_factor=row.temp_factor,
                        element=row.element,
                        charge=row.charge
                    )
                )
            pdb_lines.append(f"{'ENDMDL':<80}")
            models_str.append("\n".join(pdb_lines))
        return models_str

    def record_header(self) -> str:
        """
        PDB-formatted string representation of the record.
        """
        header = self._struct.header
        classification = "/".join(", ".join(entry) for entry in header.classification)
        dep_date = _fields.Date.to_pdb(header.dep_date)
        return f"HEADER{'':4}{classification:<40}{dep_date}{'':3}{header.pdb_id}{'':14}"

    def record_title(self) -> str:
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

    def record_compnd(self) -> str:
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

    @staticmethod
    def line_coordinates(
            is_std: bool,
            serial: int,
            atom_name: str,
            alt_loc: str,
            res_name: str,
            chain_id: str,
            res_num: int,
            res_icode: str,
            x: float,
            y: float,
            z: float,
            occupancy,
            temp_factor,
            element,
            charge
    ):
        record_name = 'ATOM  ' if is_std else 'HETATM'
        charge = '  ' if np.isnan(charge) else charge

        atom_name = f" {atom_name:<3}" if len(atom_name) < 4 else f"{atom_name:<4}"

        return (
            f"{record_name}{serial:>5} {atom_name}{alt_loc:1}"
            f"{res_name:>3} {chain_id:1}{res_num:>4}{res_icode:1}{'':3}{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{occupancy:>6.2f}{temp_factor:>6.2f}{'':10}{element.upper():>2}{charge:<2}"
        )


class EnsemblePDBWriter:

    def __init__(self, ensemble):
        self._ensemble: oc.chem.ensemble.ChemicalEnsemble = ensemble
        self._char_tab = None
        self._idx_coord_lines = None
        return

    def write(self, models, separate_models: bool = True):
        if self._char_tab is None:
            self._create_template()
        if separate_models:
            if isinstance(models, _typing.ArrayLike):
                return tuple(self._write_model(model) for model in models)
            return self._write_model(models)
        if isinstance(models, _typing.ArrayLike):
            return "\n".join(
                [self._write_model(model, end=False) for model in models[:-1]]
                + [self._write_model(models[-1])]
            )
        return self._write_model(models)

    def _write_model(self, model_num, end: bool = True):
        self._char_tab[0, 10:14] = np.char.mod("%4d", np.array([model_num + 1])).view(dtype=(str, 1))
        self._char_tab[self._idx_coord_lines, 30:54] = np.char.mod(
            '%8.3f', self._ensemble.conformation.points[model_num]
        ).view(dtype=(str, 1)).reshape(-1, 24)
        final_tab = self._char_tab if end else self._char_tab[:-1]
        return "\n".join(final_tab.view(dtype=(str, 80)).reshape(-1))

    def _create_template(self):
        atoms = self._ensemble.composition.reset_index(drop=True)
        num_polymer_chains = atoms.chain_id[atoms.res_poly].nunique()
        num_lines = atoms.shape[0] + num_polymer_chains + 3
        char_tab = np.full(shape=(num_lines, 80), fill_value=" ", dtype=(str, 1))
        char_tab[0, :5] = ["M", "O", "D", "E", "L"]
        char_tab[-2, :6] = ["E", "N", "D", "M", "D", "L"]
        char_tab[-1, :3] = ["E", "N", "D"]
        polymer_chain_ids = atoms.chain_id[atoms.res_poly]
        idx_o_terminus_atoms = polymer_chain_ids.groupby(polymer_chain_ids).aggregate(
            lambda ser: ser.last_valid_index()
        )
        idx_ter_records = idx_o_terminus_atoms + np.arange(2, idx_o_terminus_atoms.size + 2)
        o_terminus_atoms = atoms.loc[idx_o_terminus_atoms]
        char_tab[idx_ter_records, :3] = ["T", "E", "R"]
        char_tab[idx_ter_records, 6:11] = np.char.mod("%5d", o_terminus_atoms.serial + 1).view(dtype=(str, 1)).reshape(-1, 5)
        char_tab[idx_ter_records, 17:20] = o_terminus_atoms.res_name.to_numpy().astype(str).view(dtype=(str, 1)).reshape(-1, 3)
        char_tab[idx_ter_records, 21] = o_terminus_atoms.chain_id.to_numpy().astype(str)
        char_tab[idx_ter_records, 22:26] = np.char.mod("%4d", o_terminus_atoms.res_num).view(dtype=(str, 1)).reshape(-1, 4)
        char_tab[idx_ter_records, 26] = np.char.mod("%1s", o_terminus_atoms.res_icode.to_numpy(dtype=str))
        idx_coord_lines = np.setdiff1d(np.arange(1, num_lines-2), idx_ter_records)
        char_tab[idx_coord_lines, :6] = np.where(atoms.res_std, "ATOM  ", "HETATM").view(dtype=(str, 1)).reshape(-1, 6)
        char_tab[idx_coord_lines, 6:11] = np.char.mod("%5d", atoms.serial).view(dtype=(str, 1)).reshape(-1, 5)
        char_tab[idx_coord_lines, 12:16] = np.where(
            np.char.str_len(atoms.atom_name.to_numpy(dtype=str)) < 4,
            np.char.add(" ", np.char.mod("%-3s", atoms.atom_name)),
            np.char.mod("%-4s", atoms.atom_name)
        ).view(dtype=(str, 1)).reshape(-1, 4)
        char_tab[idx_coord_lines, 16] = np.char.mod("%1s", atoms.alt_loc.to_numpy(dtype=str))
        char_tab[idx_coord_lines, 17:20] = atoms.res_name.to_numpy(str).view(dtype=(str, 1)).reshape(-1, 3)
        char_tab[idx_coord_lines, 21] = atoms.chain_id.to_numpy(dtype=str)
        char_tab[idx_coord_lines, 22:26] = np.char.mod("%4d", atoms.res_num).view(dtype=(str, 1)).reshape(-1, 4)
        char_tab[idx_coord_lines, 26] = np.char.mod("%1s", atoms.res_icode.to_numpy(dtype=str))
        char_tab[idx_coord_lines, 54:60] = np.char.mod("%6.2f", atoms.occupancy).view(dtype=(str, 1)).reshape(-1, 6)
        char_tab[idx_coord_lines, 60:66] = np.char.mod("%6.2f", atoms.temp_factor).view(dtype=(str, 1)).reshape(-1, 6)
        char_tab[idx_coord_lines, 76:78] = np.char.mod("%2s", np.char.upper(atoms.element.to_numpy(dtype=str))).view(dtype=(str, 1)).reshape(-1, 2)
        char_tab[idx_coord_lines, 78:80] = np.where(
            np.isnan(atoms.charge),
            "  ",
            np.char.add(
                np.char.mod("%1.0f", np.abs(atoms.charge)),
                np.where(atoms.charge > 0, "+", "-")
            ).astype((str, 2))
        ).view(dtype=(str, 1)).reshape(-1, 2)
        self._char_tab = char_tab
        self._idx_coord_lines = idx_coord_lines
        return char_tab
