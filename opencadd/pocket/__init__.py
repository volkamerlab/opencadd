from .pocket import BindingPocket
from typing import Optional, Tuple, Sequence, Union, List, Dict
import opencadd as oc
from opencadd._http_request import HTTPRequestRetryConfig
from . import ligsite, dogsite


def from_dogsite(
        receptor,
        model: Optional[int] = 0,
        chain_id: Optional[str] = None,
        ligand_id_chain_num: Optional[Tuple[str, str, int]] = None,
        include_subpockets: bool = True,
        calculate_druggability: bool = True,
        retry_config: Optional[HTTPRequestRetryConfig] = HTTPRequestRetryConfig(),
):
    return dogsite._from_ensemble(
        ensemble=receptor,
        model=model,
        chain_id=chain_id,
        ligand_id_chain_num=ligand_id_chain_num,
        include_subpockets=include_subpockets,
        calculate_druggability=calculate_druggability,
        retry_config=retry_config
    )


# def by_ligsite(
#         receptor: oc.chem.system.ChemicalEnsemble,
#         resolution_or_grid: Union[float, Sequence[float], oc.spacetime.grid.Grid],
# ):
#     return ligsite.LigSiteDetector(receptor=receptor, resolution_or_grid=resolution_or_grid)