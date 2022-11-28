from typing import List, Union
import nglview as nv
import numpy as np


def add_points(
        view,
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
        view,
        coords: List[float],
        colors: List[float],
        radii: Union[List[float], float] = 0.2,
        opacity: float = 0.7,
):
    view._js(
        f"""
        var params = {dict(position=coords, color=colors, radius=radii)};
        var shape = new NGL.Shape('grid');
        var buffer = new NGL.SphereBuffer(params);

        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer", {{opacity:{opacity}}});
        """
    )
    return

