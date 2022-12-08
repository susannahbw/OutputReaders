from MolecularPropertyObjects import ExcitedState
from MolecularPropertyObjects import hartree_per_ev
from FileHandlingObjects import OpenFile
from MolecularPropertyObjects import VibMode
import numpy as np


class FreqOutputFile:

    def __init__(self, filepath):
        self.filepath = filepath
        self.n_atoms = None
        self.n_vib_dof = None
        self.point_group = None
        self.vib_modes = None

    def read_vib_modes(self):
        with open(self.filepath) as fp:
            file = OpenFile(fp)

            cells = file.read_to_line_containing('NAtoms').split()
            self.n_atoms = int(cells[1])

            cells = file.read_to_line_containing('Full point group').split()
            self.point_group = cells[3]

            file.read_to_line_containing('incident light, reduced masses (AMU), ')
            file.read_line()

            n_modes_read = 0
            n_dof = self.n_atoms * 3
            self.n_vib_dof = n_dof - 5 if 'inf' in self.point_group else n_dof - 6
            self.vib_modes = []
            while n_modes_read <= self.n_vib_dof:
                line = file.read_line()
                if 'Harmonic frequencies' in line:
                    break
                else:
                    cells = line.split()
                    for i, idx in enumerate(cells):
                        self.vib_modes.append(VibMode())
                        mode_ix = n_modes_read + i
                        assert int(idx) == mode_ix + 1, "Mode counting out of sync"
                        self.vib_modes[mode_ix].index = int(idx)
                        self.vib_modes[mode_ix].norm_coord_xyz = np.zeros((self.n_atoms, 3))

                    cells = file.read_line().split()
                    for i, sym in enumerate(cells):
                        mode_ix = n_modes_read + i
                        assert self.vib_modes[mode_ix].index == mode_ix + 1, \
                            "Mode counting out of sync while reading symmetries"
                        self.vib_modes[mode_ix].symmetry = sym

                    for prop in ['freq_in_inv_cm', 'red_mass_amu', 'force_const_mDyne_per_A',
                                 'IR_intensity_KM_per_mole']:
                        cells = file.read_line().split()
                        offset = [i for i, j in enumerate(cells) if j == '---'][0] + 1
                        for i, val in enumerate(cells[offset:]):
                            mode_ix = n_modes_read + i
                            assert self.vib_modes[mode_ix].index == mode_ix + 1, \
                                f"Mode counting out of sync while reading {prop}"
                            self.vib_modes[mode_ix].__setattr__(prop, float(val))

                    file.read_line()
                    for n in range(n_dof):
                        cells = file.read_line().split()
                        atom_ix = int(cells[1]) - 1
                        q = int(cells[0]) - 1
                        for i, val in enumerate(cells[3:]):
                            mode_ix = n_modes_read + i
                            assert self.vib_modes[mode_ix].index == mode_ix + 1, \
                                f"Mode counting out of sync while reading coordinates"
                            self.vib_modes[mode_ix].norm_coord_xyz[atom_ix, q] = float(val)

                    n_modes_read = len(self.vib_modes)

            for mode in self.vib_modes:
                mode.disp_vec_r = np.linalg.norm(mode.norm_coord_xyz, axis=1)

            file.f.close()


class TddftOutputFile:

    def __init__(self, filepath):
        self.filepath = filepath
        self.start_time = None
        self.end_time = None
        self.cpu_time = None
        self.input = None
        self.input_geometry = None
        self.basis = None
        self.functional = None
        self.molecular_point_group = None
        self.n_alpha_elec = None
        self.n_beta_elec = None
        self.n_basis_functions = None
        self.n_excited_states = None
        self.excited_states = None
        self.gs_energy = None

    def read_file(self):

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            file.read_to_line_containing('******************************************')
            file.read_to_line_containing('******************************************')
            line = file.read_line()
            while not ('Charge' in line) and ('Multiplicity' in line):
                self.input += line

                if 'Leave Link' in line:
                    cells = line.split()
                    self.start_time = {'day': cells[6],
                                       'month': cells[5],
                                       'year': cells[8],
                                       'time': cells[7]}

                if '#' in line:
                    cells = line.split()
                    for cell in cells:
                        if '/' in cell:
                            slash_ix = cell.index('/')
                            self.functional = cell[:slash_ix]
                            self.basis = cell[slash_ix + 1:]

                line = file.read_line()

            read_geom(file, 'input')
            self.molecular_point_group = read_quantity(file, 'Full point group', 3)
            self.n_basis_functions = int(read_quantity(file, 'basis functions', 0))
            cells = file.read_to_line_containing('alpha electrons').split()
            self.n_alpha_elec = int(cells[0])
            self.n_beta_elec = int(cells[3])

            file.read_to_line_containing('Ground to excited state transition electric dipole moments (Au):')
            file.read_line()
            line = file.read_line()
            self.excited_states = []
            while 'Ground to excited state' not in line:
                cells = line.split()
                state = ExcitedState()
                state.trans_dip = np.array([float(cells[1]),
                                            float(cells[2]),
                                            float(cells[3])])
                state.osc_strength = float(cells[5])
                self.excited_states.append(state)
                line = file.read_line()
            self.n_excited_states = len(self.excited_states)

            file.read_to_line_containing('Excitation energies and oscillator strengths:')
            for i in range(self.n_excited_states):
                cells = file.read_to_line_containing('Excited State').split()
                assert(int(cells[2][:-1]) == i + 1, 'Problem indexing excited states')
                dash_ix = cells[3].index('-')
                self.excited_states[i].spin = cells[3][:dash_ix]
                self.excited_states[i].symmetry = cells[3][dash_ix + 1:]
                self.excited_states[i].excitation_energy_eV = float(cells[4])
                self.excited_states[i].excitation_energy_au = self.excited_states[i].excitation_energy_eV \
                                                              * hartree_per_ev

                orbs_read = False
                self.excited_states[i].coeffs = []
                self.excited_states[i].from_states = []
                self.excited_states[i].to_states = []
                while not orbs_read:
                    cells = file.read_line().split()
                    if len(cells) == 0 or cells[0] == 'SavETr:':
                        orbs_read = True
                    elif cells[0] == 'This':
                        cells = file.read_to_line_containing('Total Energy').split()
                        self.gs_energy = float(cells[4]) - self.excited_states[i].excitation_energy_au
                        orbs_read = True
                    else:
                        self.excited_states[i].coeffs.append(float(cells[3]))
                        self.excited_states[i].from_states.append(int(cells[0]))
                        self.excited_states[i].to_states.append(int(cells[2]))


            read_end_of_file(file, self)


def read_geom(file: OpenFile, geom_type: str):
    """

    :param file: OpenFile object being read
    :param geom_type: Can take values 'input' or 'standard'
    :return: geom_dict: dictionary listing each atom in the structure with type and coordinates
    """

    geom_dict = {}

    assert geom_type in ['input', 'standard'], "Please choose 'input' or 'standard' as geom_type"
    if geom_type == 'input':
        search_str = 'Input orientation: '
    else:
        search_str = 'Standard orientation: '

    file.read_to_line_containing(search_str)
    for i in range(5):
        line = file.read_line()

    while not (' ---------------------------------------------------------------------'):
        cells = line.split()
        geom_type[int(cells[0])] = {'Atom type': cells[1],
                                    'x': float(cells[3]),
                                    'y': float(cells[4]),
                                    'z': float(cells[5])}
        line = file.read_line()
    return geom_dict


def read_quantity(file: OpenFile, search_str, cell_ix):
    """

    :param file: OpenFile object being read
    :param search_str: String to search for in the output file
    :param cell_ix: Index of the cell containing the desired value, once the line containing the quantity has been
    identified and split into cells. (Remember that python indexing starts from 0)
    :return: quantity (as a string)
    """

    cells = file.read_to_line_containing(search_str).split()
    return cells[cell_ix]


def read_end_of_file(file: OpenFile, results_object):
    cells = file.read_to_line_containing('Job cpu time').split()
    results_object.cpu_time = float(cells[3]) * 24 * 3600 + float(cells[5]) * 3600 \
                              + float(cells[7]) * 60 + float(cells[9])
    cells = file.read_to_line_containing('Normal termination of Gaussian')
    results_object.end_time = {'day': cells[8],
                               'month': cells[7],
                               'year': cells[10],
                               'time': cells[9]}

    file.f.close()
    return
