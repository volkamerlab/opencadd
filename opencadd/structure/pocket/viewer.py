"""
opencadd.structure.pocket.viewer

Visualize pocket(s).
"""

from matplotlib import colors
import nglview


class Viewer:
    def show(self, pocket, show_pocket_centroid=True, show_anchor_residues=True):
        """
        Visualize the pocket (subpockets, regions, and anchor residues).

        Parameters
        ----------
        show_pocket_centroid : bool
            Show the pocket centroid as sphere (default) or not.
        show_anchor_residues : bool
            Show the anchor residues as spheres (default) or not.

        Returns
        -------
        nglview.widget.NGLWidget
            Pocket visualization.
        """

        # Load structure from text in nglview
        structure = nglview.adaptor.TextStructure(pocket._text, ext=pocket._extension)
        view = nglview.widget.NGLWidget(structure)
        view._remote_call("setSize", target="Widget", args=["1000px", "600px"])
        view.clear()

        # Get residue ID > nglview index mapping
        # based on file format (this is important to know how nglview will index residues)
        residue_id2ix = self._map_residue_id2ix(pocket)

        # Show regions
        scheme_regions_list = []
        for index, region in pocket.regions.iterrows():
            color = region["region.color"]
            residue_id = region["residue.id"]
            residue_ngl_ix = residue_id2ix.loc[residue_id]
            scheme_regions_list.append([color, residue_ngl_ix])

        scheme_regions = nglview.color._ColorScheme(scheme_regions_list, label="scheme_regions")
        view.add_representation("cartoon", selection="protein", color=scheme_regions)

        # Show pocket centroids
        if show_pocket_centroid:
            view.shape.add_sphere(list(pocket.centroid), [0, 0, 1], 2, "centroid")

        # Show subpockets
        for index, subpocket in pocket.subpockets.iterrows():
            center = list(subpocket["subpocket.center"])
            name = subpocket["subpocket.name"]
            color_rgb = colors.to_rgb(subpocket["subpocket.color"])
            view.shape.add_sphere(center, color_rgb, 2, name)

        # Show anchor points
        if show_anchor_residues:
            for index, anchor_residue in pocket.anchor_residues.iterrows():
                center = list(anchor_residue["anchor_residue.center"])
                color_rgb = colors.to_rgb(anchor_residue["subpocket.color"])
                view.shape.add_sphere(center, color_rgb, 0.5)

        # Show
        return view

    def _map_residue_id2ix(self, pocket):
        """
        Map residue IDs to nglview indices depending on file format.
        In case of mol2 files, nglview will use indices starting from 1.
        In case of pdb files, nglview will use the residue IDs as indices.

        Parameters
        ----------
        file_format : string TODO
            Structural protein data format.

        Returns
        -------
        pandas.Series
            Residue IDs (index) and residue nglview indices (values).
        """

        # Get all residue names
        residue_id2ix = pocket.data[["residue.name", "residue.id"]].drop_duplicates()

        if pocket._extension == "mol2":

            # Map residue names to nglview index (starting from 1)
            residue_id2ix["residue.ngl_ix"] = [str(i) for i in range(1, len(residue_id2ix) + 1)]
            # Cast to Series (residue IDs as index, NGL index as values)
            residue_id2ix.set_index("residue.id", inplace=True)
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        else:

            # In this case, residue ID and nglview index are the same
            residue_id2ix["residue.ngl_ix"] = residue_id2ix["residue.id"]
            residue_id2ix.index = residue_id2ix["residue.id"]
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        return residue_id2ix
