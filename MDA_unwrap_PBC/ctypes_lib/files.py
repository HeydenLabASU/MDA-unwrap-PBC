"""
Location of  files
======================

Use as ::

    from MDA_unwrap_PBC.ctypes_lib.files import *

"""

__all__ = [
    "CLIB",
]

from importlib import resources
from pathlib import Path

_data_ref = resources.files("MDA_unwrap_PBC.ctypes_lib")

CLIB = (_data_ref / "unwrap.so").as_posix()

del resources
