class ExcitedState:

    def __init__(self):
        self.excitation_energy_eV = None
        self.excitation_energy_au = None
        self.symmetry = None
        self.spin = None
        self.from_states = None
        self.to_states = None
        self.coeffs = None
        self.trans_dip = None
        self.osc_strength = None


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

            line = file.read_to_line_containing('Q-Chem begins on')
            cells = line.split()
            self.start_time = {'day': cells[5],
                               'month': cells[4],
                               'year': cells[7],
                               'time': cells[6]}

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

            cells = file.read_to_line_containing('Molecular Point Group').split()
            self.molecular_point_group = cells[3]

            cells = file.read_to_line_containing('beta electrons').split()
            self.n_alpha_elec = int(cells[2])
            self.n_beta_elec = int(cells[5])

            line = file.read_to_line_containing('Requested basis set is')
            cells = line.split()
            self.basis = cells[4]
            cells = file.read_line().split()
            self.n_basis_functions = int(cells[5])

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

                            state.from_states = [cells[0] + from_spin[i][0] if spin_mask[i] and sym_mask[i] and ix_mask[i]
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

            cells = file.read_to_line_containing('Total job time:').split()
            self.cpu_time = float(cells[4][:-6])
            self.wall_time = float(cells[3][:-8])
            cells = file.read_line().split()
            self.end_time = {'day': cells[2],
                             'month': cells[1],
                             'year': cells[4],
                             'time': cells[3]}

            file.read_to_line_containing('*  Thank you very much for using Q-Chem.  Have a nice day.  *')
            file.f.close()


class OpenFile:

    def __init__(self, open_file, count=0):
        self.f = open_file
        self.count = count

    def read_to_line_containing(self, target_str):
        found = False
        while not found:
            self.count += 1
            line = self.f.readline()

            if target_str in line:
                return line
            if not line:
                break

        return 'not found'

    def read_line(self):
        self.count += 1
        return self.f.readline()
