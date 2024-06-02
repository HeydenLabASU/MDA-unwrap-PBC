"""
Test run for the MDA_unwrap_PBC package.
"""

import ctypes as ct
import numpy as np
import MDAnalysis as mda
import time

import MDA_unwrap_PBC
from MDA_unwrap_PBC.data.files import *

u = mda.Universe(TOPOLOGY,TRAJECTORY)
nRead=len(u.trajectory)

trees = unwrap.buildTrees(u)

#loop over trajectory, apply unwrap, and write out coordinates
with mda.Writer('data/traj_pbc.trr', len(u.atoms)) as w:
    for ts in u.trajectory:
        u.atoms.positions=unwrap.unwrap(u,trees)
        w.write(u)
