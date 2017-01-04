import os
import pickle
import numpy as np

from pymatgen import Structure
from abitools import AbinitTask, MrgddbTask, MrgdvTask

# =========================================================================== #
""" Tasks generation                                                       """
# =========================================================================== #


def get_den_task(dirname, structure):
    """Get a calculation for the ground state density."""

    basis_set = {
        'ecut' : 40.,
        'ecutsm' : .5,
        }
    
    kpt_grid = {
        'kptopt' : 1,
        'ngkpt' : [8, 8 , 8],
        'nshiftk' : 1,
        'shiftk' : [0.0, 0.0, 0.0],
        }
    
    
    gstate = {
        'iscf' : 5,
        'nstep' : 30,
        'tolvrs' : 1e-10,
        }

    task = AbinitTask(dirname)
    task.set_structure(structure)

    task.set_variables(basis_set)
    task.set_variables(kpt_grid)
    task.set_variables(gstate)

    return task


def get_dvscf_task_qpt(dirname, structure, qpt, den_fname, commensurate=True):
    """Get the response function task for one q-point."""

    natom = structure.num_sites

    basis_set = {
        'ecut' : 40.,
        'ecutsm' : .5,
        'nband' : 10,
        'nsym' : 1,
        }
    
    kpt_grid = {
        'kptopt' : 3,
        'ngkpt' : [8, 8, 8],
        'nshiftk' : 1,
        'shiftk' : [0.0, 0.0, 0.0],
        }
    
    gstate = {
        'iscf' : 5,
        'nstep' : 30,
        'tolvrs' : 1e-10,
        }

    dfpt = {
        'irdden' : 1,
        'irdwfk' : 1,
        'irdwfq' : 1,
        'iscf' : 7,
        'tolvrs' : 1e-16,
        'nstep' : 60,
        'rfphon' : 1,
        'rfatpol' : [1, natom],
        'rfdir' : [1, 1, 1],
        'qptopt' : 0,
        'nqpt' : 1,
        'qpt' : qpt,
        }

    wfk = {
        'irdden' : 1,
        'iscf' : -2,
        'tolwfr' : 1e-18,
        'nnsclo' : 40,
        }

    wfq = {
        'irdden' : 1,
        'iscf' : -2,
        'tolwfr' : 1e-18,
        'nnsclo' : 40,
        'nqpt' : 1,
        'qpt' : qpt,
        }

    task = AbinitTask(dirname)
    task.set_structure(structure)

    task.set_variables(basis_set)
    task.set_variables(kpt_grid)

    task.set_variables(wfk, 1)
    task.set_variables(wfq, 2)
    task.set_variables(dfpt, 3)

    if commensurate:
        task.ndtset = 2
        task.jdtset = [1,3]

        wfk_fname = task.get_odat('WFK', 1)
        wfq_fname = wfk_fname

    else:
        task.ndtset = 3
        task.jdtset = [1,2,3]

        wfk_fname = task.get_odat('WFK', 1)
        wfq_fname = task.get_odat('WFQ', 2)


    # Link files
    task.link_idat(den_fname, 'DEN', 1)
    task.link_idat(den_fname, 'DEN', 2)
    task.link_idat(den_fname, 'DEN', 3)
    task.link_idat(wfk_fname, 'WFK', 3)
    task.link_idat(wfq_fname, 'WFQ', 3)

    return task


def get_gkk_task(dirname, structure, qpt, ddb_ngqpt, den_fname,
                 ddb_fname, dvdb_fname, commensurate=False):
    """Get the electron-phonon coupling task."""

    natom = structure.num_sites

    basis_set = {
        'ecut' : 40.,
        'ecutsm' : .5,
        'nband' : 14,
        }
    
    kpt_grid = {
        'kptopt' : 3,
        'ngkpt' : [8, 8, 8],
        'nshiftk' : 1,
        'shiftk' : [0.0, 0.0, 0.0],
        }

    gkk = {
        'optdriver' : 7,
        'eph_task' : 2,
        'prtphdos' : 0,
        'irdden' : 1,
        'irdwfk' : 1,
        'irdwfq' : 1,
        'iscf' : -2,
        'tolwfr' : 1e-18,
        'qptopt' : 0,
        'nqpt' : 1,
        'qpt' : qpt,
        'ddb_ngqpt' : ddb_ngqpt,
        'ddb_shiftq' : 3*[.0],
        }

    wfk = {
        'irdden' : 1,
        'iscf' : -2,
        'tolwfr' : 1e-18,
        'nnsclo' : 40,
        'istwfk' : '*1',
        }

    wfq = {
        'irdden' : 1,
        'iscf' : -2,
        'tolwfr' : 1e-18,
        'nnsclo' : 40,
        'nqpt' : 1,
        'qpt' : qpt,
        'istwfk' : '*1',
        }

    task = AbinitTask(dirname)
    task.set_structure(structure)

    task.set_variables(basis_set)
    task.set_variables(kpt_grid)

    task.set_variables(wfk, 1)
    task.set_variables(wfq, 2)
    task.set_variables(gkk, 3)

    if commensurate:
        task.ndtset = 2
        task.jdtset = [1,3]

        wfk_fname = task.get_odat('WFK', 1)
        wfq_fname = wfk_fname

    else:
        task.ndtset = 3
        task.jdtset = [1,2,3]

        wfk_fname = task.get_odat('WFK', 1)
        wfq_fname = task.get_odat('WFQ', 2)


    # Link files
    task.link_idat(den_fname, 'DEN', 1)
    task.link_idat(den_fname, 'DEN', 2)
    task.link_idat(den_fname, 'DEN', 3)

    task.link_idat(wfk_fname, 'WFK', 3)
    task.link_idat(wfq_fname, 'WFQ', 3)

    task.link_idat(ddb_fname, 'DDB', 3)
    task.link_idat(dvdb_fname, 'DVDB', 3)

    return task


def get_mrgddb_task(dirname, dvscf_tasks):
    """Make a task for merging the DDB."""
    ddb_fnames = [task.get_odat('DDB', 3) for task in dvscf_tasks]
    return MrgddbTask(dirname, ddb_fnames)


def get_mrgdv_task(dirname, dvscf_tasks, natom):
    """Make a task for merging the POT files."""
    pot_files = list()
    for i, task in enumerate(dvscf_tasks):
        for iat in range(natom):
            for icart in range(3):
                ipert = 3*iat + icart + 1
                pot_files.append(task.get_odat('POT{}'.format(ipert), 3))

    return MrgdvTask(dirname, pot_files)


# =========================================================================== #
""" Execution                                                               """
# =========================================================================== #


# Open the Structure from a file
structure_fname = '../Data/Structures/LiF/data-1-1-LiF-relaxed.json'
structure = Structure.from_file(structure_fname)
natom = structure.num_sites


# Get the list of q-points from a file.
qptgrid_fname = '../Data/Kptgrids/LiF/Symmetrized/kpt-g4.pkl'
with open(qptgrid_fname, 'r') as f:
    qptgrid = pickle.load(f)
    ngqpt = qptgrid['ngkpt']
    qpts = qptgrid['kpt']
    nqpt = len(qpts)


# List of q-points onto which the potential will be interpolated
qpts_gkk = [[0.0125, 0.0, 0.0]]


# Set up all tasks
topdir = 'Gkk-workflow'

# Density calculation
den_task = get_den_task(os.path.join(topdir, 'Density'), structure)
den_fname = den_task.get_odat('DEN')


# DVSCF calculations
dvscf_tasks = list()
for iqpt, qpt in enumerate(qpts):

    subdir = 'DVSCF-qpt-{:0=3}'.format(iqpt+1)
    dirname = os.path.join(topdir, subdir)
    
    task = get_dvscf_task_qpt(dirname, structure, qpt, den_fname)
    
    dvscf_tasks.append(task)


# mrgddb task
mrgddb_task = get_mrgddb_task(os.path.join(topdir, 'Mrgddb'), dvscf_tasks)
merged_ddb_fname = mrgddb_task.ddb_fname


# mrgdv task
mrgdv_task = get_mrgdv_task(os.path.join(topdir, 'Mrgdv'), dvscf_tasks, natom)
merged_dvdb_fname = mrgdv_task.dvdb_fname


# GKK tasks
gkk_tasks = list()
for iqpt, qpt in enumerate(qpts_gkk):

    subdir = 'GKK-qpt-{:0=3}'.format(iqpt+1)
    dirname = os.path.join(topdir, subdir)

    task = get_gkk_task(dirname, structure, qpt, ngqpt, den_fname,
                        merged_ddb_fname, merged_dvdb_fname)
    
    gkk_tasks.append(task)


# Set binaries and 
all_tasks = [den_task] + dvscf_tasks + [mrgddb_task, mrgdv_task] + gkk_tasks
for task in all_tasks:

    task.set_bindir('/path/to/abinit/build/src/98_main/')

    task.pseudo_dir = '../Data/Pseudos'
    task.pseudos = ['03-Li.LDA.TM.pspnc', '09-F.LDA.TM.pspnc']


# Set some variables
for task in ([den_task] + dvscf_tasks + gkk_tasks):

    task.nproc = 4
    #task.nodes = 1
    #task.nproc_per_node = 32

    task.set_variables({
        'ecut' : 10.,
        'ngkpt' : [4, 4, 4],
        })

for task in gkk_tasks:

    task.set_variables({
        'ngkpt' : [2, 2, 2],
        'nband' : 14,
        })


# Write files and report
for task in all_tasks:
    task.write()
    #task.run()
    task.report()

