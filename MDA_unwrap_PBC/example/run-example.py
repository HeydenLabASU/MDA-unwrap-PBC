"""
Usage example for mda_unwrap_pbc package
"""

import MDAnalysis as mda
import MDA_unwrap_PBC as pbc
from MDA_unwrap_PBC.example.datafiles import *

#input parameters

u = mda.Universe(TOPOL,TRAJ)
nRead=len(u.trajectory)

unwrap = pbc.unwrap(u)

#loop over trajectory, apply unwrap, and write out coordinates
with mda.Writer(TRAJ_PBC, len(u.atoms)) as w:
    for ts in u.trajectory:
        unwrap.single_frame()
        w.write(u)
