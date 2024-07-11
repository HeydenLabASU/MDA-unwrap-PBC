import ctypes as ct
import numpy as np
import MDAnalysis as mda
import time
from MDA_unwrap_PBC.ctypes_lib.files import *

class t_branch(ct.Structure):
    """recursive datatype to describe molecules as trees"""
    pass
t_branch._fields_ = (("node", ct.c_int32),
                   ("links",ct.POINTER(t_branch)),
                   ("nLinks",ct.c_int32))

class t_trees(ct.Structure):
    """array of molecule trees"""
    _fields_ = (("trees",ct.POINTER(t_branch)),
                 ("nTrees",ct.c_int32))

#load the shared library with C routines
clib = ct.cdll.LoadLibrary(CLIB)

#define argument types of function 'buildTrees' in imported library 'clib'
clib.buildTrees.argtypes = [
    ct.POINTER(t_trees),
    ct.c_int32,
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ct.c_int32,
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
]
#define return type of function 'buildTrees' in imported library 'clib'
clib.buildTrees.restype = ct.c_int32

#define argument types of function 'unwrap' in imported library 'clib'
clib.unwrap.argtypes = [
    t_trees,
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
]
#define return type of function 'unwrap' in imported library 'clib'
clib.unwrap.restype = ct.c_int32

class unwrap:
    """functions to make molecules whole in PBC trajectories"""
    def buildTrees(u):
        """build recursive bond trees that define molecules"""
        start=time.process_time()
        nAtoms=len(u.atoms)
        atomTags=np.zeros(nAtoms,dtype=np.int32)

        bondList=np.zeros([len(u.bonds),2],dtype=np.int32)
        i=0
        for b in u.bonds:
            bondList[i][0]=b[0].index
            bondList[i][1]=b[1].index
            i+=1
        nBonds=i
        bondTags=np.zeros(len(u.bonds),dtype=np.int32)

        masses=u.atoms.masses.astype(np.float32)

        trees = t_trees()
        error = clib.buildTrees(
            ct.pointer(trees),
            ct.c_int(nAtoms),
            atomTags,
            ct.c_int(nBonds),
            bondTags,
            bondList,
            masses
        )
        if error != 0:
            print(f'ERROR reported by \'buildTrees\' function\nsee \'error.log\'\n')
        stop=time.process_time()
        print(f'detected {trees.nTrees} molecules')
        print(f'tree building time: {stop-start:.2f}s')
        return trees

    def unwrap(u,trees):
        """
        unwrap coordinates: 
        ensure that all components of covalent bonds are shorter than 1/2 the box
        """
        coords=u.atoms.positions.astype(np.float32)
        error = clib.unwrap(
            trees,
            coords,
            u.dimensions
        )
        if error != 0:
            print(f'ERROR reported by \'unwrap\' function\nsee \'error.log\'\n')
        return coords
