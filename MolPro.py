from MolecularPropertyObjects import CasExcitedState
from MolecularPropertyObjects import CasPlusExcitedState
from FileHandlingObjects import OpenFile
import numpy as np


class CasOutputFile:

    def __init__(self, filepath, extension=None):
        self.filepath = filepath
        self.start_time = None
        self.total_calc_wall_time_s = None
        self.cas_cpu_time = None
        self.memory_mb = None
        self.input = None
        self.geometry = None
        self.basis = None
        self.molecular_point_group = None
        self.n_excited_states = None
        self.gs_energy_au = None
        self.excited_states = None
        self.orbs = None
        if extension:
            assert extension in ['caspt2'], "Extension must be one of ['caspt2']"
        self.extension = extension
        if self.extension == 'caspt2':
            self.rs2_cpu_time = None

    def read_file(self):
        with open(self.filepath) as fp:
            file = OpenFile(fp)

            file.read_to_line_containing('***,')
            self.input = ''
            input_end = False
            while not input_end:
                line = file.read_line()
                if 'Variables initialized' in line:
                    input_end = True
                elif not line:
                    break
                else:
                    self.input += line

            read_start_time(file, self)
            self.basis = read_property(file, 'SETTING BASIS', 3)
            self.molecular_point_group = read_property(file, 'Point group', 2)
            read_structure(file, self)

            file.read_to_line_containing('1PROGRAM * MULTI')

            line = file.read_to_line_containing('State symmetry')
            sym_dict = {}
            state_symmetries_read = False
            while not state_symmetries_read:
                cells = line.split()
                ix = int(cells[2])
                sym_dict[ix] = {}
                cells = file.read_to_line_containing('Number of electrons:').split()
                sym_dict[ix]['nelec'] = int(cells[3])
                sym_dict[ix]['spin'] = cells[5][9:]
                sym_dict[ix]['space sym'] = cells[7][9:]
                cells = file.read_to_line_containing('Number of states:').split()
                sym_dict[ix]['nstates'] = int(cells[3])
                file.read_line()
                file.read_line()
                line = file.read_line()
                if 'State symmetry' not in line:
                    state_symmetries_read = True

            self.n_excited_states = 0
            self.excited_states = []
            for ix in sym_dict.keys():
                for j in range(sym_dict[ix]['nstates']):
                    cells = file.read_to_line_containing('Results for state').split()
                    if self.extension in ['caspt2']:
                        state = CasPlusExcitedState()
                    else:
                        state = CasExcitedState()
                    self.n_excited_states += 1
                    state.symmetry = cells[3][(cells[3].index('.') + 1):]
                    assert state.symmetry == sym_dict[ix]['space sym'], 'Problem reading state symmetry'
                    state.spin = sym_dict[ix]['spin']

                    cells2 = file.read_to_line_containing(f'!MCSCF STATE {cells[3]} Energy').split()
                    if cells2[2] == '1.1':
                        self.gs_energy_au = {'CAS': float(cells2[4])}
                        state.excited_state_energy_au = float(cells2[4])
                    else:
                        state.excited_state_energy_au = float(cells2[4])

                    state.config = []
                    state.coeffs = []

                    self.excited_states.append(state)

            state_ix = 0
            for ix in sym_dict.keys():
                file.read_to_line_containing('CI vector for state symmetry')
                file.read_line()
                file.read_line()
                ci_vec_read = False
                while not ci_vec_read:
                    line = file.read_line()
                    if line == '\n':
                        ci_vec_read = True
                    else:
                        cells = line.split()
                        for j in range(sym_dict[ix]['nstates']):
                            self.excited_states[state_ix + j].config.append(cells[0] + ' ' + cells[1])
                            self.excited_states[state_ix + j].coeffs.append(float(cells[2 + j]))
                state_ix += sym_dict[ix]['nstates']

            if self.extension == 'caspt2':
                read_rs2_output(file, self, sym_dict)

            cells = file.read_to_line_containing('PROGRAMS   *').split()
            cas_ix = [i for i in range(len(cells)) if cells[i] == 'CASSCF'][0] + 1
            if self.extension == 'caspt2':
                rs_ix = [i for i in range(len(cells)) if 'RS2' in cells[i]][0] + 1
            cells = file.read_line().split()
            self.cas_cpu_time = float(cells[cas_ix])
            if self.extension == 'caspt2':
                self.rs2_cpu_time = float(cells[rs_ix])
            cells = file.read_line().split()
            self.total_calc_wall_time_s = float(cells[3])
            cells = file.read_line().split()
            self.memory_mb = float(cells[3])


            line = file.read_to_line_containing('DUMP ORBITAL')
            orbs_read = False
            self.orbs = []
            while not orbs_read:
                cells = line.split()
                orb = Orbital()
                orb.index = int(cells[5])
                orb.sym_ix = cells[2]
                orb.occupation = float(cells[7])
                orb.eigenvalue = float(cells[9])
                self.orbs.append(orb)

                line = file.read_line()
                if line == '\n':
                    orbs_read = True

            file.f.close()


def read_start_time(file: OpenFile, results_object):
    line = file.read_to_line_containing('DATE')
    cells = line.split()
    results_object.start_time = {'day': cells[5][:2],
                                 'month': cells[5][3:6],
                                 'year': cells[5][-2:],
                                 'time': cells[7]}


def read_property(file: OpenFile, search_str: str, cell_ix_for_quantity: int):
    cells = file.read_to_line_containing(search_str).split()
    return cells[cell_ix_for_quantity]


def read_structure(file: OpenFile, results_object):
    file.read_to_line_containing('ATOMIC COORDINATES')
    for i in range(3):
        file.read_line()

    results_object.geometry = []
    geometry_read = False
    while not geometry_read:
        line = file.read_line()
        if line == '\n':
            geometry_read = True
        else:
            cells = line.split()
            results_object.geometry.append({'atom' : cells[1],
                                            'xyz'  : np.array([float(cells[3]),
                                                               float(cells[4]),
                                                               float(cells[5])])})
    return


def read_rs2_output(file: OpenFile, results_object, sym_dict):
    file.read_to_line_containing('1PROGRAM * RS2')

    cells = file.read_to_line_containing('Level shift=').split()
    for i in range(results_object.n_excited_states):
        results_object.excited_states[i].level_shift = float(cells[2])

    i = 0
    for ix in sym_dict.keys():
        cells = file.read_to_line_containing('RESULTS FOR STATE').split()
        sym = int(cells[3][:cells[3].index('.')])
        if sym != ix:
            continue
        for j in range(sym_dict[ix]['nstates']):
            if str(ix) + '.' + str(j + 1) == cells[3]:
                cells = file.read_to_line_containing('!RSPT2 STATE').split()
                if cells[2] == '1.1':
                    results_object.gs_energy_au['CASPT2'] = float(cells[4])

                results_object.excited_states[i].cas_excited_state_energy_au = \
                    results_object.excited_states[i].excited_state_energy_au * 1
                results_object.excited_states[i].excited_state_energy_au = float(cells[4])
                i += 1
            else:
                i += 1
                continue
    return

class Orbital:

    def __init__(self):
        self.index = None
        self.sym_ix = None
        self.occupation = None
        self.eigenvalue = None
