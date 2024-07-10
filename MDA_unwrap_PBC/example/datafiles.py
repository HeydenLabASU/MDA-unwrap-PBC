"""
Location of data files
======================

Use as ::

    from MDA_unwrap_PBC.example.datafiles import *

"""

__all__ = [
    "TOPOL",
    "TRAJ",
    "TRAJ_PBC",
]

from importlib import resources
from pathlib import Path

_data_ref = resources.files("MDA_unwrap_PBC.data")

TOPOL = (_data_ref / "topol.tpr").as_posix()
TRAJ = (_data_ref / "traj.trr").as_posix()
TRAJ_PBC = (_data_ref / "traj_pbc.trr").as_posix()

del resources
