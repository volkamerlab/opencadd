from typing import List, Union, NoReturn, Optional,Sequence
import nglview as nv
import numpy as np
from opencadd._typing import ArrayLike


class NGLViewAdaptor(nv.Structure, nv.Trajectory):

    def __init__(self, ensemble):
        self._ensemble = ensemble
        self.ext = "pdb"
        self.params = {}
        self.id = 0
        return

    def get_structure_string(self):
        return self._ensemble.to_pdb(0)

    def get_coordinates(self, index):
        return self._ensemble.conformation.points[index]

    @property
    def n_frames(self):
        return self._ensemble.conformation.points.shape[0]


class NGLViewer:

    def __init__(self):
        self._widget = nv.NGLWidget()
        return

    @property
    def nglwidget(self):
        return self._widget

    def add_ensemble(self, ensemble):
        self._widget.add_trajectory(NGLViewAdaptor(ensemble=ensemble))
        return

    def add_points(
            self,
            coords: List[float],
            colors: List[float],
            opacity: float = 0.7,
    ):
        self._widget._js(
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
            self,
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
            raise ValueError("Size of `radii` and `coords` do not match.")

        self._widget._js(
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
            self,
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
        self._widget._execute_js_code(
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
            self,
            component_id: int = 0,
            selection: str = "protein",
            representation: str = "cartoon",
            aspect_ratio: float = 1,
            scale: float = 1,
            multiple_bond: bool = True,
    ):
        self._widget._execute_js_code(
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

    def remove_component_by_name(self, name: str) -> NoReturn:
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
        self._widget._execute_js_code(
            f"""this.stage.removeComponent(this.stage.getComponentsByName("{name}").first)"""
        )
        return
