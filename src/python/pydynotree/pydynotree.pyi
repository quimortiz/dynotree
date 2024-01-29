"""

        pydynotree
        -----------------------

        .. currentmodule:: pydynotree


"""

from __future__ import annotations
import numpy

__all__ = [
    "R2",
    "R2SO2",
    "R3",
    "R3SO3",
    "R4",
    "R5",
    "R6",
    "R7",
    "RX",
    "SO2",
    "SO3",
    "SpaceX",
    "TreeR2",
    "TreeR2SO2",
    "TreeR2SO2dp",
    "TreeR2dp",
    "TreeR3SO3",
    "TreeR3SO3dp",
    "TreeR4",
    "TreeR4dp",
    "TreeR7",
    "TreeR7dp",
    "TreeRX",
    "TreeRXdp",
    "TreeSO2",
    "TreeSO2dp",
    "TreeSO3",
    "TreeSO3dp",
    "TreeX",
    "TreeXdp",
    "rand",
    "rand01",
    "srand",
]

class R2:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[2, 1]],
        arg1: numpy.ndarray[numpy.float64[2, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[2, 1]],
        arg1: numpy.ndarray[numpy.float64[2, 1]],
        arg2: numpy.ndarray[numpy.float64[2, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[2, 1]],
        arg1: numpy.ndarray[numpy.float64[2, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[2, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[2, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[2, 1]],
        arg1: numpy.ndarray[numpy.float64[2, 1]],
    ) -> None: ...

class R2SO2:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
        arg2: numpy.ndarray[numpy.float64[3, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[2, 1]],
        arg1: numpy.ndarray[numpy.float64[2, 1]],
    ) -> None: ...

class R3:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
        arg2: numpy.ndarray[numpy.float64[3, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
    ) -> None: ...

class R3SO3:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
        arg2: numpy.ndarray[numpy.float64[7, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[7, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[7, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[3, 1]],
        arg1: numpy.ndarray[numpy.float64[3, 1]],
    ) -> None: ...

class R4:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
        arg2: numpy.ndarray[numpy.float64[4, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
    ) -> None: ...

class R5:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[5, 1]],
        arg1: numpy.ndarray[numpy.float64[5, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[5, 1]],
        arg1: numpy.ndarray[numpy.float64[5, 1]],
        arg2: numpy.ndarray[numpy.float64[5, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[5, 1]],
        arg1: numpy.ndarray[numpy.float64[5, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[5, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[5, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[5, 1]],
        arg1: numpy.ndarray[numpy.float64[5, 1]],
    ) -> None: ...

class R6:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[6, 1]],
        arg1: numpy.ndarray[numpy.float64[6, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[6, 1]],
        arg1: numpy.ndarray[numpy.float64[6, 1]],
        arg2: numpy.ndarray[numpy.float64[6, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[6, 1]],
        arg1: numpy.ndarray[numpy.float64[6, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[6, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[6, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[6, 1]],
        arg1: numpy.ndarray[numpy.float64[6, 1]],
    ) -> None: ...

class R7:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
        arg2: numpy.ndarray[numpy.float64[7, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[7, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[7, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[7, 1]],
        arg1: numpy.ndarray[numpy.float64[7, 1]],
    ) -> None: ...

class RX:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
        arg2: numpy.ndarray[numpy.float64[m, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[m, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[m, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
    ) -> None: ...

class SO2:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[1, 1]],
        arg1: numpy.ndarray[numpy.float64[1, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[1, 1]],
        arg1: numpy.ndarray[numpy.float64[1, 1]],
        arg2: numpy.ndarray[numpy.float64[1, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[1, 1]],
        arg1: numpy.ndarray[numpy.float64[1, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[1, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[1, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[1, 1]],
        arg1: numpy.ndarray[numpy.float64[1, 1]],
    ) -> None: ...

class SO3:
    def __init__(self) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
        arg2: numpy.ndarray[numpy.float64[4, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: numpy.ndarray[numpy.float64[4, 1]],
        arg1: numpy.ndarray[numpy.float64[4, 1]],
    ) -> None: ...

class SpaceX:
    def __init__(self, arg0: list[str]) -> None: ...
    def distance(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
    ) -> float: ...
    def distance_to_rectangle(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
        arg2: numpy.ndarray[numpy.float64[m, 1]],
    ) -> float: ...
    def interpolate(
        self,
        arg0: numpy.ndarray[numpy.float64[m, 1]],
        arg1: numpy.ndarray[numpy.float64[m, 1]],
        arg2: float,
        arg3: numpy.ndarray[numpy.float64[m, 1], numpy.ndarray.flags.writeable],
    ) -> None: ...
    def sample_uniform(
        self, arg0: numpy.ndarray[numpy.float64[m, 1], numpy.ndarray.flags.writeable]
    ) -> None: ...
    def set_bounds(
        self,
        arg0: list[numpy.ndarray[numpy.float64[m, 1]]],
        arg1: list[numpy.ndarray[numpy.float64[m, 1]]],
    ) -> None: ...

class TreeR2:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[2, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> R2: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: R2 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[2, 1]]) -> TreeR2dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[2, 1]], arg1: float
    ) -> list[TreeR2dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[2, 1]], arg1: int
    ) -> list[TreeR2dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeR2SO2:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> R2SO2: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: R2SO2 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[3, 1]]) -> TreeR2SO2dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: float
    ) -> list[TreeR2SO2dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: int
    ) -> list[TreeR2SO2dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeR2SO2dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeR2dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeR3SO3:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> R3SO3: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: R3SO3 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[7, 1]]) -> TreeR3SO3dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: float
    ) -> list[TreeR3SO3dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: int
    ) -> list[TreeR3SO3dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeR3SO3dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeR4:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> R4: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: R4 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[4, 1]]) -> TreeR4dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: float
    ) -> list[TreeR4dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: int
    ) -> list[TreeR4dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeR4dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeR7:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> R7: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: R7 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[7, 1]]) -> TreeR7dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: float
    ) -> list[TreeR7dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[7, 1]], arg1: int
    ) -> list[TreeR7dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeR7dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeRX:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> RX: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: RX = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[m, 1]]) -> TreeRXdp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: float
    ) -> list[TreeRXdp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: int
    ) -> list[TreeRXdp]: ...
    def splitOutstanding(self) -> None: ...

class TreeRXdp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeSO2:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[1, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> SO2: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: SO2 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[1, 1]]) -> TreeSO2dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[1, 1]], arg1: float
    ) -> list[TreeSO2dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[1, 1]], arg1: int
    ) -> list[TreeSO2dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeSO2dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeSO3:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> SO3: ...
    def init_tree(
        self, runtime_dimension: int = -1, t_state_space: SO3 = ...
    ) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[4, 1]]) -> TreeSO3dp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: float
    ) -> list[TreeSO3dp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[4, 1]], arg1: int
    ) -> list[TreeSO3dp]: ...
    def splitOutstanding(self) -> None: ...

class TreeSO3dp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

class TreeX:
    def __init__(self) -> None: ...
    def addPoint(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: int, arg2: bool
    ) -> None: ...
    def getStateSpace(self) -> SpaceX: ...
    def init_tree(self, arg0: int, arg1: SpaceX) -> None: ...
    def search(self, arg0: numpy.ndarray[numpy.float64[m, 1]]) -> TreeXdp: ...
    def searchBall(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: float
    ) -> list[TreeXdp]: ...
    def searchKnn(
        self, arg0: numpy.ndarray[numpy.float64[m, 1]], arg1: int
    ) -> list[TreeXdp]: ...
    def splitOutstanding(self) -> None: ...

class TreeXdp:
    def __init__(self) -> None: ...
    @property
    def distance(self) -> float: ...
    @property
    def id(self) -> int: ...

def rand() -> int: ...
def rand01() -> float: ...
def srand(arg0: int) -> None: ...
