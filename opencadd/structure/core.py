"""
Core objects used across the library
"""

from MDAnalysis import Universe
import numpy as np


class Structure(Universe):
    """
    Core object to load structures with.

    Thin wrapper around MDAnalysis.Universe objects
    """

    @classmethod
    def from_pdbid(cls, pdbid):
        import mmtf

        return cls(mmtf.fetch(pdbid))

    @classmethod
    def from_string(cls, pdbid_or_path, **kwargs):
        import os

        if os.path.isfile(pdbid_or_path):
            return cls(pdbid_or_path, **kwargs)
        return cls.from_pdbid(pdbid_or_path, **kwargs)

    @classmethod
    def from_atomgroup(cls, *args):
        """Create a new new :class:`Structure` from one or more
        ``AtomGroup`` instances.

        Parameters
        ----------
        *args : ``AtomGroup``
            One or more AtomGroups.

        Returns
        -------
        structure : :class:`Structure`

        Raises
        ------
        ValueError
            Too few arguments or an AtomGroup is empty and
        TypeError
            Arguments are not :class:`AtomGroup` instances.

        Notes
        -----
        This is take straight from ``MDAnalysis.universe``. Refer to that
        module for more information.

        """
        from MDAnalysis.coordinates.memory import MemoryReader
        from MDAnalysis.topology.base import squash_by
        from MDAnalysis.core import groups
        from MDAnalysis.core.topology import Topology
        from MDAnalysis.core.topologyattrs import AtomAttr, ResidueAttr, SegmentAttr

        if len(args) == 0:
            raise ValueError("Need at least one AtomGroup for merging")

        for a in args:
            if not isinstance(a, groups.AtomGroup):
                raise TypeError(repr(a) + " is not an AtomGroup")
        for a in args:
            if len(a) == 0:
                raise ValueError("cannot merge empty AtomGroup")

        # Create a new topology using the intersection of topology attributes
        blank_topology_attrs = set(dir(Topology(attrs=[])))
        common_attrs = set.intersection(*[set(dir(ag.universe._topology)) for ag in args])
        tops = set(["bonds", "angles", "dihedrals", "impropers"])

        attrs = []

        # Create set of attributes which are array-valued and can be simply
        # concatenated together
        common_array_attrs = common_attrs - blank_topology_attrs - tops
        # Build up array-valued topology attributes including only attributes
        # that all arguments' universes have
        for attrname in common_array_attrs:
            for ag in args:
                attr = getattr(ag, attrname)
                attr_class = type(getattr(ag.universe._topology, attrname))
                if issubclass(attr_class, AtomAttr):
                    pass
                elif issubclass(attr_class, ResidueAttr):
                    attr = getattr(ag.residues, attrname)
                elif issubclass(attr_class, SegmentAttr):
                    attr = getattr(ag.segments, attrname)
                else:
                    raise NotImplementedError(
                        "Don't know how to handle"
                        " TopologyAttr not subclassed"
                        " from AtomAttr, ResidueAttr,"
                        " or SegmentAttr."
                    )
                if type(attr) != np.ndarray:
                    raise TypeError(
                        "Encountered unexpected topology "
                        "attribute of type {}".format(type(attr))
                    )
                try:
                    attr_array.extend(attr)
                except NameError:
                    attr_array = list(attr)
            attrs.append(attr_class(np.array(attr_array, dtype=attr.dtype)))
            del attr_array

        # Build up topology groups including only those that all arguments'
        # universes have
        for t in tops & common_attrs:
            offset = 0
            bondidx = []
            types = []
            for ag in args:
                # create a mapping scheme for this atomgroup
                mapping = {a.index: i for i, a in enumerate(ag, start=offset)}
                offset += len(ag)

                tg = getattr(ag, t)
                bonds_class = type(getattr(ag.universe._topology, t))
                # Create a topology group of only bonds that are within this ag
                # ie we don't want bonds that extend out of the atomgroup
                tg = tg.atomgroup_intersection(ag, strict=True)

                # Map them so they refer to our new indices
                new_idx = [tuple([mapping[x] for x in entry]) for entry in tg.indices]
                bondidx.extend(new_idx)
                if hasattr(tg, "_bondtypes"):
                    types.extend(tg._bondtypes)
                else:
                    types.extend([None] * len(tg))
            if any(t is None for t in types):
                attrs.append(bonds_class(bondidx))
            else:
                types = np.array(types, dtype="|S8")
                attrs.append(bonds_class(bondidx, types))

        # Renumber residue and segment indices
        n_atoms = sum([len(ag) for ag in args])
        residx = []
        segidx = []
        res_offset = 0
        seg_offset = 0
        for ag in args:
            # create a mapping scheme for this atomgroup's parents
            res_mapping = {r.resindex: i for i, r in enumerate(ag.residues, start=res_offset)}
            seg_mapping = {r.segindex: i for i, r in enumerate(ag.segments, start=seg_offset)}
            res_offset += len(ag.residues)
            seg_offset += len(ag.segments)

            # Map them so they refer to our new indices
            residx.extend([res_mapping[x] for x in ag.resindices])
            segidx.extend([seg_mapping[x] for x in ag.segindices])

        residx = np.array(residx, dtype=np.int32)
        segidx = np.array(segidx, dtype=np.int32)

        _, _, [segidx] = squash_by(residx, segidx)

        n_residues = len(set(residx))
        n_segments = len(set(segidx))

        top = Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )

        # Create and populate a universe
        coords = np.vstack([a.positions for a in args])
        u = cls(top, coords[None, :, :], format=MemoryReader)

        return u

    def write(self, *args, **kwargs):
        # Workaround for https://github.com/MDAnalysis/mdanalysis/issues/2679
        if self.dimensions is None:
            self.trajectory.ts._unitcell = np.zeros(6)
        return self.atoms.write(*args, **kwargs)
