from typing import Literal, Optional, Union, Sequence, Tuple
import jax.numpy as jnp
import pandas as pd

import opencadd as oc
import opencadd.io.pdb._writer


class ChemicalEnsemble:

    def __init__(self, composition: pd.DataFrame, conformation: jnp.ndarray):
        self._composition = composition
        self._conformation = oc.spacetime.pointcloud.DynamicPointCloud(data=conformation)
        self._pdb_writer = oc.io.pdb._writer.EnsemblePDBWriter(ensemble=self)
        return

    @property
    def composition(self):
        return self._composition

    @property
    def conformation(self):
        return self._conformation

    def remove(self, *args: Literal["nonpoly"]):
        composition = self._composition[self._composition.res_poly]
        conformation = self._conformation.points[:, self._composition.res_poly.to_numpy()]
        return ChemicalEnsemble(composition=composition, conformation=conformation)

    def to_pdb(
            self, model: Union[int, Sequence[int]],
            separate_models: bool = True
    ) -> Union[str, Tuple[str, ...]]:
        return self._pdb_writer.write(models=model, separate_models=separate_models)


def from_pdb_structure(structure):
    atoms = structure.atom
    composition = atoms#.drop(["model_num", "alt_loc", "occupancy", "x", "y", "z", "temp_factor"], axis=1)
    conformation = jnp.expand_dims(jnp.array(atoms[["x", "y", "z"]]), axis=0)
    return ChemicalEnsemble(composition=composition, conformation=conformation)


def from_pdb_id(
        pdb_id: str,
        biological_assembly_id: Optional[int] = 1,
):
    return from_pdb_structure(
        oc.io.pdb.read.from_pdb_id(pdb_id=pdb_id, biological_assembly_id=biological_assembly_id)
    )
