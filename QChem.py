from MolecularPropertyObjects import ExcitedState
from FileHandlingObjects import OpenFile


class EomEeCcsdOutputFile:

    def __init__(self, filepath):
        self.filepath = filepath
        self.start_time = None
        self.end_time = None
        self.wall_time = None
        self.cpu_time = None
        self.input = None
        self.geometry = None
        self.basis = None
        self.molecular_point_group = None
        self.n_alpha_elec = None
        self.n_beta_elec = None
        self.n_basis_functions = None
        self.n_excited_states = None
        self.excited_states = None

    def read_file(self):

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            read_start_time(file, self)

            file.read_to_line_containing('User input:')
            file.read_line()

            input_end = False
            self.input = ''
            state_sets = []
            while not input_end:
                line = file.read_line()
                if '------------' in line:
                    input_end = True
                elif not line:
                    break
                else:
                    self.input += line
                    if 'EE_SINGLETS' in line or 'EE_TRIPLETS' in line:
                        cells = line.split()
                        cells = cells[2:]
                        cells[0] = cells[0][1:]
                        cells[-1] = cells[-1][:-1]
                        for cell in cells[:-1]:
                            state_sets.append(int(cell[:-1]))
                        state_sets.append(int(cells[-1]))
            self.n_excited_states = sum(state_sets)

            read_point_group(file, self)
            read_n_electrons(file, self)
            read_n_basis_functions(file, self)

            self.excited_states = []
            for n in range(len(state_sets)):
                cells = file.read_to_line_containing('Solving for EOMEE-CCSD').split()
                symm = cells[3]
                spin = cells[4]

                for m in range(state_sets[n]):
                    state = ExcitedState()
                    state.symmetry = symm
                    state.spin = spin

                    cells = file.read_to_line_containing('Excitation energy').split()
                    state.excitation_energy_eV = float(cells[8])

                    file.read_to_line_containing('Transitions between orbitals')
                    orbs_read = False
                    state.coeffs = []
                    from_ix = []
                    from_sym = []
                    from_spin = []
                    to_ix = []
                    to_sym = []
                    to_spin = []
                    while not orbs_read:
                        cells = file.read_line().split()
                        if len(cells) == 0:
                            orbs_read = True
                        else:
                            state.coeffs.append(float(cells[0]))
                            from_ix.append(int(cells[1]))
                            from_sym.append(cells[2])
                            from_spin.append('Alpha' if (cells[3] == 'A') else 'Beta')
                            to_ix.append(int(cells[5]))
                            to_sym.append(cells[6])
                            to_spin.append('Alpha' if (cells[7] == 'A') else 'Beta')

                    state.from_states = [0 for i in range(len(from_ix))]
                    state.to_states = [0 for i in range(len(from_ix))]
                    file.read_to_line_containing('Number')
                    orbs_read = False
                    while not orbs_read:
                        cells = file.read_line().split()
                        if len(cells) == 0:
                            break

                        if cells[1] == 'Occ':
                            spin_mask = [(i == cells[2]) for i in from_spin]
                            sym_mask = [(i == cells[4]) for i in from_sym]
                            ix_mask = [(i == int(cells[3])) for i in from_ix]

                            state.from_states = [
                                cells[0] + from_spin[i][0] if spin_mask[i] and sym_mask[i] and ix_mask[i]
                                else state.from_states[i] for i in range(len(state.from_states))]
                        elif cells[1] == 'Vir':
                            spin_mask = [(i == cells[2]) for i in to_spin]
                            sym_mask = [(i == cells[4]) for i in to_sym]
                            ix_mask = [(i == int(cells[3])) for i in to_ix]

                            state.to_states = [cells[0] + to_spin[i][0] if spin_mask[i] and sym_mask[i] and ix_mask[i]
                                               else state.to_states[i] for i in range(len(state.to_states))]
                        else:
                            print('error reading state indices')
                            orbs_read = True

                    self.excited_states.append(state)

            read_end_of_file(file, self)


class TddftOutputFile:

    def __init__(self, filepath):
        self.filepath = filepath
        self.start_time = None
        self.end_time = None
        self.wall_time = None
        self.cpu_time = None
        self.input = None
        self.geometry = None
        self.basis = None
        self.molecular_point_group = None
        self.n_alpha_elec = None
        self.n_beta_elec = None
        self.n_basis_functions = None
        self.n_excited_states = None
        self.excited_states = None

    def read_file(self, read_tda=False):

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            read_start_time(file, self)

            file.read_to_line_containing('User input:')
            file.read_line()

            input_end = False
            self.input = ''
            n = 0
            while not input_end:
                line = file.read_line()
                if '------------' in line:
                    input_end = True
                elif not line:
                    break
                else:
                    self.input += line
                    if 'CIS_N_ROOTS' in line:
                        cells = line.split()
                        nroots = int(cells[2])

                    if 'CIS_SINGLETS' in line or 'CIS_TRIPLETS' in line:
                        cells = line.split()
                        n = n + 1 if cells[2] == 'true' else n

                    if 'RPA' in line:
                        rpa = True if 'true' in line else False
                        read_tda *= rpa

            self.n_excited_states = n * nroots

            read_point_group(file, self)
            read_n_electrons(file, self)
            read_n_basis_functions(file, self)

            if read_tda:
                file.read_to_line_containing('TDDFT/TDA Excitation Energies')
                shift = 0
            else:
                file.read_to_line_containing('TDDFT Excitation Energies')
                shift = 1

            self.excited_states = []
            for m in range(self.n_excited_states):
                state = ExcitedState()
                cells = file.read_to_line_containing('Excited state').split()
                state.excitation_energy_eV = float(cells[7])

                cells = file.read_to_line_containing('Multiplicity').split()
                state.spin = cells[1]

                cells = file.read_to_line_containing('Trans. Mom').split()
                state.trans_dip = [float(cells[i]) for i in [2, 4, 6]]

                cells = file.read_to_line_containing('Strength').split()
                state.osc_strength = float(cells[2])

                orbs_read = False
                state.coeffs = []
                state.from_states = []
                state.to_states = []
                while not orbs_read:
                    cells = file.read_line().split()
                    if len(cells) == 0:
                        orbs_read = True
                    else:
                        state.coeffs.append(float(cells[shift + 7]))
                        state.from_states.append(int(cells[shift + 1][:-1]))
                        state.to_states.append(int(cells[shift + 4][:-1]) + self.n_alpha_elec)

                self.excited_states.append(state)

            read_end_of_file(file, self)


def read_start_time(file: OpenFile, results_object):
    line = file.read_to_line_containing('Q-Chem begins on')
    cells = line.split()
    results_object.start_time = {'day': cells[5],
                                 'month': cells[4],
                                 'year': cells[7],
                                 'time': cells[6]}
    return


def read_point_group(file: OpenFile, results_object):
    cells = file.read_to_line_containing('Molecular Point Group').split()
    results_object.molecular_point_group = cells[3]
    return


def read_n_electrons(file: OpenFile, results_object):
    cells = file.read_to_line_containing('beta electrons').split()
    results_object.n_alpha_elec = int(cells[2])
    results_object.n_beta_elec = int(cells[5])
    return


def read_n_basis_functions(file: OpenFile, results_object):
    line = file.read_to_line_containing('Requested basis set is')
    cells = line.split()
    results_object.basis = cells[4]
    cells = file.read_line().split()
    results_object.n_basis_functions = int(cells[5])
    return


def read_end_of_file(file: OpenFile, results_object):
    cells = file.read_to_line_containing('Total job time:').split()
    results_object.cpu_time = float(cells[4][:-6])
    results_object.wall_time = float(cells[3][:-8])
    cells = file.read_line().split()
    results_object.end_time = {'day': cells[2],
                               'month': cells[1],
                               'year': cells[4],
                               'time': cells[3]}

    file.read_to_line_containing('*  Thank you very much for using Q-Chem.  Have a nice day.  *')
    file.f.close()
    return
