"""
Usage example for mda_unwrap_pbc package
"""

import MDAnalysis as mda
import MDA_unwrap_PBC as unwrap
from MDA_unwrap_PBC.example.datafiles import *

#input parameters

u = mda.Universe(TOPOL,TRAJ)
nRead=len(u.trajectory)

trees = unwrap.unwrap.buildTrees(u)

#loop over trajectory, apply unwrap, and write out coordinates
with mda.Writer(TRAJ_PBC, len(u.atoms)) as w:
    for ts in u.trajectory:
        u.atoms.positions=unwrap.unwrap.unwrap(u,trees)
        w.write(u)
