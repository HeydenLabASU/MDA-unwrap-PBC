"""
Location of data files
======================

Use as ::

    from MDA_unwrap_PBC.data.files import *

"""

__all__ = [
    "TOPOLOGY",  # GROMACS topology file for small protein in water
    "TRAJECTORY",   # GROMACS trajectory with PBC applied
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

from pkg_resources import resource_filename

TOPOLOGY = resource_filename(__name__, "topol.tpr")
TRAJECTORY  = resource_filename(__name__, "traj.trr")
MDANALYSIS_LOGO = resource_filename(__name__, "mda.txt")
