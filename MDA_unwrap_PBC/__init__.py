"""
MDA-unwrap-PBC
Fast C routines to make molecules whole in PBC trajectories
"""

# Add imports here
from .unwrap import *
from importlib.metadata import version

__version__ = version("MDA-unwrap-PBC")
