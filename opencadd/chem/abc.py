from abc import ABC, abstractmethod
from typing import Sequence, Optional, Literal

from opencadd._typing import ArrayLike


class Composition(ABC):

    @abstractmethod
    def __init__(
            self,
            atomic_numbers: Sequence[int],
            mass_numbers: Optional[Sequence[int]] = None,
            charges: Optional[Sequence[float]] = None,
    ):
        ...

    @abstractmethod
    @property
    def atomic_numbers(self) -> ArrayLike:
        ...

    @abstractmethod
    @property
    def atomic_masses(self) -> ArrayLike:
        ...

    @abstractmethod
    @property
    def charges(self) -> ArrayLike:
        ...

    @abstractmethod
    def van_der_waals_radii(self, ref):
        ...


class Configuration(ABC):
    ...


class Conformation(ABC):
    ...


class Molecule(ABC):

    @abstractmethod
    @property
    def composition(self) -> Composition:
        ...

    @abstractmethod
    @property
    def configuration(self) -> Configuration:
        ...

    @abstractmethod
    @property
    def conformation(self) -> Conformation:
        ...


class Polymer(Molecule, ABC):
    ...


class Protein(Polymer, ABC):
    pass


class ProteinComposition(Composition):

    @abstractmethod
    @property
    def atom_types_autodock(self) -> ArrayLike:
        ...

    @abstractmethod
    @property
    def residue_type_per_atom(self) -> ArrayLike:
        ...