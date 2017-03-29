"""
The users define their own workflow.
This can be done through class definition or through simple scripts.
"""

# =========================================================================== #
""" Main function                                                           """
# =========================================================================== #

from pymatgen import Structure
import pickle

def main():

    flow = ElectronPhononSelfEnergyFlow(

        dirname = 'Calculation',
        structure = Structure.from_file('../Data/Structures/LiF/data-1-1-LiF-relaxed.json'),
        pseudo_dir = '../Data/Pseudos',
        pseudos = ['03-Li.LDA.TM.pspnc', '09-F.LDA.TM.pspnc'],

        qpts = [
            [ 0.00000000E+00, 0.00000000E+00, 0.00000000E+00],
            [ 5.00000000E-01, 0.00000000E+00, 0.00000000E+00],
            [ 5.00000000E-01, 5.00000000E-01, 0.00000000E+00],
            ],

        variables = {
            'ecut' : 5.0,
            'nband' : 8,
            'ixc' : 1,
            'enunit' : 2,
            'diemac' : 9.0,
            'npfft' : 1,
            },

        den_variables = {
            'ngkpt' : [4, 4, 4],
            'tolvrs' : 1e-16,
            },

        wfk_dvscf_variables = {
            'ngkpt' : [4, 4, 4],
            'tolwfr' : 1e-14,
            },

        wfk_epc_variables = {
            'ngkpt' : [1, 1, 1],
            'tolwfr' : 1e-14,
            },

        dvscf_variables = {
            'ngkpt' : [4, 4, 4],
            'tolvrs' : 1e-10,
            },

        epc_variables = {
            'ngkpt' : [1, 1, 1],
            'tolwfr' : 1e-14,
            },

        )

    flow.write()
    #flow.run()
    flow.report()


# =========================================================================== #
""" Class definition                                                       """
# =========================================================================== #

import os
from abitools import Workflow, AbinitTask, MrgddbTask, MrgdvTask  # 4.2

class ElectronPhononSelfEnergyFlow(Workflow):

    default_parameters = {

        # Basis set
        'ecut' : 10,

        # kpt_grid
        'kptopt' : 1,
        'ngkpt' : [2, 2 , 2],
        'nshiftk' : 1,
        'shiftk' : [0.0, 0.0, 0.0],
        }

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------

        variables: dict
            Variables used for all calculations
        den_variables: dict
            Variables used for density
        wfk_dvscf_variables: dict
            Variables used for wavefunctions for dvscf
        wfk_epc_variables: dict
            Variables used for wavefunctions for epc
        dvscf_variables: dict
            Variables used for dvscf
        epc_variables: dict
            Variables used for epc
    
        """

        super(ElectronPhononSelfEnergyFlow, self).__init__(**kwargs)

        kwargs.pop('dirname')

        self.structure = kwargs.pop('structure')

        self.natom = self.structure.num_sites
        self.ncart = 3

        self.qpts =  kwargs.get('qpts', [3*[.0]])
        self.nqpt = len(self.qpts)

        self.setup_tasks(**kwargs)

    def setup_tasks(self, **kwargs):
        """Setup all tasks for the workflow."""

        self.clear_tasks()

        # Common variables for all tasks
        variables = kwargs.pop('variables', {})

        # ------------------------------------------------------------------- #

        # Density calculation

        # Input variables
        self.den_variables = kwargs.get('den_variables', {})
        for k, d in variables.items():
            self.den_variables.setdefault(k, d)

        # Initialize task
        self.den_task = self.get_den_task(
            dirname = os.path.join(self.dirname, 'Den'),
            structure = self.structure,
            variables = self.den_variables,
            **kwargs)

        self.add_task(self.den_task)

        # ------------------------------------------------------------------- #

        # Wavefunction calculation for response function

        # Input variables
        self.wfk_dvscf_variables = kwargs.get('wfk_dvscf_variables', {})
        for k, d in variables.items():
            self.wfk_dvscf_variables.setdefault(k, d)

        # Initialize task
        self.wfk_dvscf_task = self.get_wfk_task(
            dirname = os.path.join(self.dirname, 'DVSCF', 'WFK'),
            structure = self.structure,
            den_fname = self.den_fname,
            variables = self.wfk_dvscf_variables,
            **kwargs)

        self.add_task(self.wfk_dvscf_task)

        # ------------------------------------------------------------------- #

        # Wavefunction calculation for epc

        # Input variables
        self.wfk_epc_variables = kwargs.get('wfk_epc_variables', {})
        for k, d in variables.items():
            self.wfk_epc_variables.setdefault(k, d)

        # Initialize tasks
        self.wfk_epc_task = self.get_wfk_task(
            dirname = os.path.join(self.dirname, 'EPC', 'WFK'),
            structure = self.structure,
            den_fname = self.den_fname,
            variables = self.wfk_epc_variables,
            **kwargs)

        self.add_task(self.wfk_epc_task)

        #---------------------------------------------------------------------#

        # DVSCF calculations
    
        # Input variables
        self.dvscf_variables = kwargs.get('dvscf_variables', {})
        for k, d in variables.items():
            self.dvscf_variables.setdefault(k, d)

        # Initialize tasks
        self.wfq_dvscf_tasks = list()
        self.dvscf_tasks = list()
        for iqpt, qpt in enumerate(self.qpts):

            wfq_dvscf_task = self.get_wfq_task(
                dirname = os.path.join(
                    self.dirname,
                    'DVSCF',
                    'qpt-{:0=4}'.format(iqpt+1),
                    'WFQ',
                    ),
                structure = self.structure,
                qpt = qpt,
                den_fname = self.den_fname,
                variables = self.wfk_dvscf_variables,
                **kwargs)

            self.wfq_dvscf_tasks.append(wfq_dvscf_task)
            self.add_task(wfq_dvscf_task)

            dvscf_task = self.get_dvscf_task_qpt(
                dirname = os.path.join(
                    self.dirname,
                    'DVSCF',
                    'qpt-{:0=4}'.format(iqpt+1),
                    'DVSCF',
                    ),
                structure = self.structure,
                qpt = qpt,
                den_fname = self.den_fname,
                wfk_fname = self.wfk_dvscf_fname,
                wfq_fname = wfq_dvscf_task.get_odat('WFQ'),
                variables = self.dvscf_variables,
                **kwargs)


            self.dvscf_tasks.append(dvscf_task)
            self.add_task(dvscf_task)


        #---------------------------------------------------------------------#

        # Electron-phonon coupling calculations
    
        # Input variables
        self.epc_variables = kwargs.get('epc_variables', {})
        for k, d in variables.items():
            self.epc_variables.setdefault(k, d)

        # Initialize tasks
        self.wfq_epc_tasks = list()
        self.epc_tasks = list()
        for iqpt, qpt in enumerate(self.qpts):

            wfq_epc_task = self.get_wfq_task(
                dirname = os.path.join(
                    self.dirname,
                    'EPC',
                    'qpt-{:0=4}'.format(iqpt+1),
                    'WFQ',
                    ),
                structure = self.structure,
                qpt = qpt,
                den_fname = self.den_fname,
                variables = self.wfk_epc_variables,
                **kwargs)

            self.wfq_epc_tasks.append(wfq_epc_task)
            self.add_task(wfq_epc_task)

            epc_task = self.get_epc_task(
                dirname = os.path.join(
                    self.dirname,
                    'EPC',
                    'qpt-{:0=4}'.format(iqpt+1),
                    'EPC',
                    ),
                structure = self.structure,
                qpt = qpt,
                den_fname = self.den_fname,
                wfk_fname = self.wfk_epc_fname,
                wfq_fname = wfq_epc_task.get_odat('WFQ'),
                den1_fnames = self.get_qpt_den1_fnames(iqpt),
                variables = self.epc_variables,
                **kwargs)

            self.epc_tasks.append(epc_task)
            self.add_task(epc_task)

        #---------------------------------------------------------------------#

        return

    @property
    def den_fname(self):
        return self.den_task.get_odat('DEN')

    @property
    def wfk_epc_fname(self):
        return self.wfk_epc_task.get_odat('WFK')

    @property
    def wfk_dvscf_fname(self):
        return self.wfk_dvscf_task.get_odat('WFK')

    @property
    def den1_fnames(self):
        """Return the perturbed DEN files for dvscf MassLauncher."""
        fnames = list()
        for iqpt in range(self.nqpt):
            fnames.extend(self.get_qpt_den1_fnames(iqpt))
        return fnames

    def get_qpt_den1_fnames(self, iqpt):
        """Return the perturbed DEN files for dvscf MassLauncher."""
        fnames = list()
        for iatom in range(self.natom):
            for icart in range(self.ncart):
                itask = iqpt #* natom * ncart + iatom * ncart + icart
                ipert = iatom * self.ncart + icart + 1
                task = self.dvscf_tasks[itask]
                fnames.append(task.get_odat('DEN{}'.format(ipert)))
    
        return fnames


    def get_den_task(self, dirname, structure, **kwargs):
        """Get a calculation for the ground state density."""

        variables = kwargs.get('variables', {})
    
        local_default_parameters = {
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)

        gstate = {
            'iscf' : 7,
            'nstep' : 30,
            'tolvrs' : 1e-10,
            }
    
        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(local_default_parameters)
        task.set_variables(gstate)
        task.set_variables(variables)

        return task

    def get_wfk_task(self, dirname, structure, den_fname, **kwargs):
        """Get a calculation for the wavefunctions."""
    
        wfk = {
            'irdden' : 1,
            'iscf' : -2,
            }
    
        local_default_parameters = {
            'kptopt' : 3,
            'tolwfr' : 1e-18,
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)

        variables = kwargs.get('variables', {})

        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(local_default_parameters)
        task.set_variables(wfk)
        task.set_variables(variables)

        task.link_idat(den_fname, 'DEN')

        return task    

    def get_wfq_task(self, dirname, structure, qpt, den_fname, **kwargs):
        """Get a calculation for the q-shifted wavefunctions."""

        variables = kwargs.pop('variables', {})

        variables.update({
            'nqpt' : 1,
            'qpt' : qpt,
            })

        kwargs['variables'] = variables

        return self.get_wfk_task(dirname, structure, den_fname, **kwargs)

    def get_dvscf_task_qpt(self, dirname, structure, qpt,
                           den_fname, wfk_fname, wfq_fname,
                           commensurate=False, **kwargs):
        """Get the response function task for one q-point."""
    
        local_default_parameters = {
            'kptopt' : 3,
            'nsym' : 1,
            'nband' : 10,
            'nnsclo' : 40,
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)

        variables = kwargs.get('variables', {})

        natom = structure.num_sites
        ncart = 3
    
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
    
        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(local_default_parameters)
    
        task.set_variables(dfpt)
    
        task.set_variables(variables)
    
        # Link files
        task.link_idat(den_fname, 'DEN')
        task.link_idat(wfk_fname, 'WFK')
        task.link_idat(wfq_fname, 'WFQ')

        return task

    def get_epc_task(self, dirname, structure, qpt,
                     den_fname, wfk_fname, wfq_fname, den1_fnames,
                     commensurate=False, **kwargs):
        """Compute the EIGR2D and GKK elements."""

        natom = structure.num_sites
        ncart = 3
    
        local_default_parameters = {
            'kptopt' : 3,
            'nsym' : 1,
            'nnsclo' : 40,
            'tolwfr' : 1e-18,
            'enunit' : 2,
            'optforces' : 1,
            'istwfk' : '*1',
            'autoparal' : 1,
            'kptopt' : 3,
            }

        for k, d in self.default_parameters.items():
            local_default_parameters.setdefault(k, d)

        variables = kwargs.get('variables', {})

        epc = {
            'rfphon' : 1,
            'rfatpol' : [1, natom],
            'rfdir' : [1, 1, 1],
            'iscf' : -2,
            'irdwfk' : 1,
            'irdwfq' : 1,
            'ird1den' : 1,
            'ieig2rf' : 5,   # print _GKK.nc file alongside the _EIGR2D.nc file.
            'bdeigrf' : -1,  # maximum number of bands for ZPR
            'elph2_imagden' : '0.05 eV',
            #'nsym' : 1,
            'qptopt' : 0,
            'nqpt' : 1,
            'qpt' : qpt,
            }

        task = AbinitTask(dirname, structure=structure, **kwargs)
    
        task.set_variables(local_default_parameters)
    
        task.set_variables(epc)

        task.set_variables(variables)
    
        task.link_idat(den_fname, 'DEN')
        task.link_idat(wfk_fname, 'WFK')
        task.link_idat(wfq_fname, 'WFQ')

        # Link the perturbed density previously calculated.
        for iatom in range(natom):
            for icart in range(ncart):
                ipert = iatom * ncart + icart
                task.link_idat(den1_fnames[ipert], 'DEN{}'.format(ipert+1))

        return task


# =========================================================================== #
""" Execution                                                               """
# =========================================================================== #

if __name__ == '__main__':
    main()
