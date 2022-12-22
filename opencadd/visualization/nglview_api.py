from typing import List, Union, NoReturn, Optional,Sequence
import nglview as nv
import numpy as np
from opencadd.typing import ArrayLike


def add_points(
        view: nv.widget.NGLWidget,
        coords: List[float],
        colors: List[float],
        opacity: float = 0.7,
):
    view._js(
        f"""
        var point_buffer = new NGL.PointBuffer(
            {{
                position: new Float32Array({coords}),
                color: new Float32Array({colors}),
            }},
            {{
                sizeAttenuation: true,
                pointSize: 2,
                opacity: 0.1,
                useTexture: true,
                alphaTest: 0.0,
                edgeBleach: 0.7,
                forceTransparent: true,
                sortParticles: true
            }}
        )
        var shape = new NGL.Shape('grid');

        shape.addBuffer(point_buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{ opacity: {opacity} }});
        """
    )
    return

def add_spheres(
        view: nv.widget.NGLWidget,
        coords: Sequence[float],
        colors: Sequence[float] = (0.5, 0.5, 0.5),
        radii: Union[Sequence[float], float] = 0.2,
        opacity: float = 0.7,
        name: Optional[str] = "spheres"
):
    if not isinstance(coords, np.ndarray):
        coords = np.array(coords, dtype=np.single)
    if not isinstance(colors, np.ndarray):
        colors = np.array(colors, dtype=np.half)
    if not isinstance(radii, np.ndarray):
        radii = np.array(radii, dtype=np.half)

    num_coords = coords.size
    if num_coords % 3 != 0:
        raise ValueError(
            "There must be 3n values in `coords`, corresponding to n (x,y,z)-coordinates."
        )
    num_points = num_coords // 3
    num_rgbs = colors.size
    if num_rgbs < 3:
        raise ValueError("`colors` must at least be a single RGB value.")
    elif num_rgbs == 3:
        colors = np.tile(colors, reps=num_points)
    elif num_rgbs != num_coords:
        raise ValueError("Size of `colors` and `coords` do not match.")
    num_radii = radii.size
    if num_radii == 1:
        radii = np.tile(radii, reps=num_points)
    elif num_radii != num_points:
        raise ValueError("Size of `radii` and `coords` do not mtch.")

    view._js(
        f"""
        var params = {
            dict(
                position=coords.flatten().tolist(), 
                color=colors.flatten().tolist(), 
                radius=radii.flatten().tolist()
            )
        };
        var shape = new NGL.Shape('{name}');
        var buffer = new NGL.SphereBuffer(params);

        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{opacity:{opacity}}});
        """
    )
    return


def add_representation_within_radius_of_selection(
        view: nv.widget.NGLWidget,
        component_id: int = 0,
        selection: str = "ligand",
        radius: float = 4.0,
        representation_type: str = "licorice"
) -> NoReturn:
    """
    Add a representation to all residues that have atoms within a certain radius of an atom in a
    selection.

    Parameters
    ----------
    view : nglview.widget.NGLWidget
    component_id : int
    selection : str
    radius : float
    representation_type : str

    Returns
    -------
    None
    """
    view._execute_js_code(
        f"""
        var component = this.stage.compList[{component_id}];
        var selection = new NGL.Selection("{selection}");
        var atomsWithinRad = component.structure.getAtomSetWithinSelection(selection, {radius});
        // Expand selection to all atoms within residues of selected atoms.
        var residuesWithinRad = component.structure.getAtomSetWithinGroup(atomsWithinRad);
        component.addRepresentation(
            "{representation_type}", 
            {{sele: residuesWithinRad.toSeleString()}}
        );
        """
    )

def add_representation_to_structure(
        view: nv.widget.NGLWidget,
        component_id: int = 0,
        selection: str = "protein",
        representation: str = "cartoon",
        aspect_ratio: float = 1,
        scale: float = 1,
        multiple_bond: bool = True,
):
    view._execute_js_code(
        f"""
        // Get the component
        var component = this.stage.compList[{component_id}];
        // exit if the component does not contain a structure
        if(component.type !== "structure") return;
        // add representation
        component.addRepresentation(
            "{representation}", 
            {{
                sele: {selection}, 
                aspectRatio: {aspect_ratio}, 
                scale: {scale}, 
                multipleBond: {multiple_bond}
            }} 
        );
        """
    )


def remove_component_by_name(view: nv.widget.NGLWidget, name: str) -> NoReturn:
    """
    Given a component's name, remove it from the stage of a given NGLWidget.

    Parameters
    ----------
    view : nglview.widget.NGLWidget
        The widget to remove the component from.
    name : str
        Name of the component to remove.

    Returns
    -------
    None
    """
    view._execute_js_code(
        f"""this.stage.removeComponent(this.stage.getComponentsByName("{name}").first)"""
    )
    return
