"""
Handlers used to report results of analysis
"""


def atomium_to_nglview(*models):
    """
    Represent Atomium models in NGLView

    Parameters
    ----------
    models : atomium.Model

    Returns
    -------
    nglview.NGLWidget
    """
    import nglview as nv
    from structuralalignment.utils import enter_temp_directory

    v = nv.NGLWidget()
    with enter_temp_directory():
        for i, model in enumerate(models):
            fn = f"tmp{i}.pdb"
            model.save(fn)
            fs = nv.FileStructure(fn)
            v.add_component(fs)
    return v
