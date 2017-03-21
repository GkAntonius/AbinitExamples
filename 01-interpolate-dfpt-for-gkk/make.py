
# =========================================================================== #
""" Main function                                                           """
# =========================================================================== #

from pymatgen import Structure
import pickle

def main():

    # Get the list of q-points from a file.
    qptgrid_fname = '../Data/Kptgrids/LiF/Symmetrized/kpt-g4.pkl'
    with open(qptgrid_fname, 'r') as f:
        qptgrid = pickle.load(f)
    
        ngqpt = qptgrid['ngkpt']
        qpts = qptgrid['kpt']
    
    # Initialize the workflow 
    flow = GkkFlow(
    
        dirname = 'Gkk-workflow',
    
        # Open structure file
        structure = Structure.from_file(
            '../Data/Structures/LiF/data-1-1-LiF-relaxed.json'),
    
        # Q-point grid
        ngqpt = ngqpt,
        qpts = qpts,
        qpts_gkk = [[.0125,.0,.0]],
        pseudo_dir = '../Data/Pseudos',
        pseudos = ['03-Li.LDA.TM.pspnc', '09-F.LDA.TM.pspnc'],
    
        # Specific variables
        density_extra_variables = {
            'ecut' : 10.,
            'ngkpt' : [4, 4, 4],
            },
    
        dvscf_extra_variables = {
            'ecut' : 10.,
            'ngkpt' : [4, 4, 4],
            },
    
        gkk_extra_variables = {
            'ecut' : 10.,
            'ngkpt' : [2, 2, 2],
            'nband' : 14,
            },
    
        commensurate=False,
    
        nproc = 4
    
        # Specify executable if not in environment
        #bindir = '/path/to/abinit/build/src/98_main/',
        )
    
    # Write files and report
    flow.write()
    #flow.run()
    flow.report()

# =========================================================================== #
""" Class definition                                                       """
# =========================================================================== #

import os
from abitools import Workflow, AbinitTask, MrgddbTask, MrgdvTask  # 4.2

class GkkFlow(Workflow):

    default_parameters = {

        # Basis set
        'ecut' : 10,
        'ecutsm' : .5,

        # kpt_grid
        'kptopt' : 1,
        'ngkpt' : [8, 8 , 8],
        'nshiftk' : 1,
        'shiftk' : [0.0, 0.0, 0.0],

        }

    def __init__(self, **kwargs):

        super(GkkFlow, self).__init__(**kwargs)

        kwargs.pop('dirname')

        self.structure = kwargs.pop('structure')

        self.natom = self.structure.num_sites

        self.ngqpt =  kwargs.get('ngqpt', [1,1,1])
        self.qpts =  kwargs.get('qpts', [3*[.0]])
        self.qpts_gkk =  kwargs.get('qpts_gkk', self.qpts)
        self.nqpt = len(self.qpts)

        self.setup_tasks(**kwargs)

    def setup_tasks(self, **kwargs):

        # Density calculation
        den_dirname = os.path.join(self.dirname, 'Density')
        self.den_task = self.get_den_task(
            den_dirname,
            self.structure,
            extra_variables = kwargs.get('density_extra_variables', {}),
            **kwargs)

        self.den_fname = self.den_task.get_odat('DEN')


        # DVSCF calculations
        self.dvscf_tasks = list()
        for iqpt, qpt in enumerate(self.qpts):

            dvscf_dirname = os.path.join(
                self.dirname,
                'DVSCF-qpt-{:0=3}'.format(iqpt+1)
                )

            task = self.get_dvscf_task_qpt(
                dvscf_dirname,
                self.structure,
                qpt,
                self.den_fname,
                extra_variables = kwargs.get('dvscf_extra_variables', {}),
                **kwargs)
            
            self.dvscf_tasks.append(task)

        # mrgddb task
        self.mrgddb_task = self.get_mrgddb_task(
            os.path.join(self.dirname, 'Mrgddb'),
            self.dvscf_tasks,
            **kwargs)

        self.merged_ddb_fname = self.mrgddb_task.ddb_fname

        # mrgdv task
        self.mrgdv_task = self.get_mrgdv_task(
            os.path.join(self.dirname, 'Mrgdv'),
            self.dvscf_tasks,
            self.natom,
            **kwargs)

        self.merged_dvdb_fname = self.mrgdv_task.dvdb_fname

        # GKK tasks
        self.gkk_tasks = list()
        for iqpt, qpt in enumerate(self.qpts_gkk):
        
            
            gkk_dirname = os.path.join(
                self.dirname,
                'GKK-qpt-{:0=3}'.format(iqpt+1),
                )
        
            task = self.get_gkk_task(
                gkk_dirname,
                self.structure,
                qpt,
                self.ngqpt,
                self.den_fname,
                self.merged_ddb_fname,
                self.merged_dvdb_fname,
                extra_variables = kwargs.get('gkk_extra_variables', {}),
                **kwargs)
            
            self.gkk_tasks.append(task)

        for task in (
            [self.den_task]
            + self.dvscf_tasks
            + [self.mrgddb_task, self.mrgdv_task]
            + self.gkk_tasks
            ):

            self.add_task(task)

    def get_den_task(self, dirname, structure, **kwargs):
        """Get a calculation for the ground state density."""

        extra_variables = kwargs.get('extra_variables', {})
    
        gstate = {
            'iscf' : 5,
            'nstep' : 30,
            'tolvrs' : 1e-10,
            }
    
        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(self.default_parameters)
        task.set_variables(gstate)
    
        task.set_variables(extra_variables)

        return task
    
    
    def get_dvscf_task_qpt(self, dirname, structure, qpt, den_fname,
                           commensurate=True, **kwargs):
        """Get the response function task for one q-point."""
    
        local_default_parameters = {
            'kptopt' : 3,
            'nsym' : 1,
            'nband' : 10,
            'nnsclo' : 40,
            'tolwfr' : 1e-18,
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)

        extra_variables = kwargs.get('extra_variables', {})

        natom = structure.num_sites
    
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
            }
    
        wfq = {
            'irdden' : 1,
            'iscf' : -2,
            'nqpt' : 1,
            'qpt' : qpt,
            }
    
        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(local_default_parameters)
    
        task.set_variables(wfk, 1)
        task.set_variables(wfq, 2)
        task.set_variables(dfpt, 3)
    
        task.set_variables(extra_variables)
    
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
    
    
    def get_gkk_task(self, dirname, structure, qpt, ddb_ngqpt, den_fname,
                     ddb_fname, dvdb_fname, commensurate=False, **kwargs):
        """Get the electron-phonon coupling task."""
    
        extra_variables = kwargs.get('extra_variables', {})

        local_default_parameters = {
            'kptopt' : 3,
            #'nsym' : 1,
            'nband' : 10,
            'nnsclo' : 40,
            'tolwfr' : 1e-18,
            'istwfk' : '*1',
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)


        natom = structure.num_sites
    
        gkk = {
            'optdriver' : 7,
            'eph_task' : 2,
            'prtphdos' : 0,
            'irdden' : 1,
            'irdwfk' : 1,
            'irdwfq' : 1,
            'iscf' : -2,
            'qptopt' : 0,
            'nqpt' : 1,
            'qpt' : qpt,
            'ddb_ngqpt' : ddb_ngqpt,
            'ddb_shiftq' : 3*[.0],
            }
    
        wfk = {
            'irdden' : 1,
            'iscf' : -2,
            }
    
        wfq = {
            'irdden' : 1,
            'iscf' : -2,
            'nqpt' : 1,
            'qpt' : qpt,
            }
    
        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(wfk, 1)
        task.set_variables(wfq, 2)
        task.set_variables(gkk, 3)
    
        task.set_variables(local_default_parameters)
        task.set_variables(extra_variables)

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

    def get_mrgddb_task(self, dirname, dvscf_tasks, **kwargs):
        """Make a task for merging the DDB."""
        ddb_fnames = [task.get_odat('DDB', 3) for task in dvscf_tasks]
        return MrgddbTask(dirname, ddb_fnames)


    def get_mrgdv_task(self, dirname, dvscf_tasks, natom, **kwargs):
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

if __name__ == '__main__':
    main()
