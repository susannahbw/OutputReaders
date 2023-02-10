from MolecularPropertyObjects import CasExcitedState
from MolecularPropertyObjects import CasPt2ExcitedState
from MolecularPropertyObjects import ExcitedState
from MolecularPropertyObjects import SsSrCasPt2ExcitedState
from MolecularPropertyObjects import hartree_per_ev
from FileHandlingObjects import OpenFile
from FileHandlingObjects import find_next_nondigit_character
from Gaussian import get_len_of_current_result_list
import numpy as np
import math

spin_dict = {0: 'Singlet',
             1: 'Doublet',
             2: 'Triplet',
             'Singlet': 0,
             'Doublet': 1,
             'Triplet': 2}


class MolProOutputFile:
    def __init__(self, filepath):
        self.filepath = filepath
        self.results = {}

    def read_file(self):

        readfuncs_dict = {"1PROGRAM * MULTI (Direct Multiconfiguration SCF)": read_cas_result,
                          "1PROGRAM * RS2C (Multireference RS Perturbation Theory)": read_rs2_result}
        flags = [k for k in readfuncs_dict.keys()]

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            while True:
                line = file.read_line()
                if not line:
                    break
                str_check = [f in line for f in flags]
                if any(str_check):
                    flag = [flags[i] for i in range(len(flags)) if str_check[i]]
                    assert len(flag) == 1, 'Expected only one flag per line'
                    flag = flag[0]

                    print(f"Found \'{flag}\' on line {file.count}")
                    readfuncs_dict[flag](self, file)

            file.f.close()
        return


class MolProInput:

    def __init__(self, input_str):
        self.input_str = input_str
        self.programs = None

    def split_into_programs(self):
        open_brackets = [pos for pos, char in enumerate(self.input_str) if char == '{']
        close_brackets = [pos for pos, char in enumerate(self.input_str) if char == '}']

        assert len(open_brackets) == len(close_brackets), "Found more { than } in Molpro input"
        n_progs = len(open_brackets)
        self.programs = {}

        separators = [',', ';', '\n']
        for i in range(n_progs):
            self.programs[i] = {}
            self.programs[i]['cmd_str'] = self.input_str[open_brackets[i + 1]:close_brackets[i]]
            sep_ix = [pos for pos, char in enumerate(self.programs[i]['cmd_str']) if char in separators]
            first_word = self.programs[i]['cmd_str'][:sep_ix[0]]
            self.programs[i]['type'] = first_word

        return


class CasOutputFile:

    def __init__(self, filepath, es_type=CasExcitedState):
        self.filepath = filepath
        self.start_time = None
        self.total_calc_wall_time_s = None
        self.cas_cpu_time = None
        self.memory_mb = None
        self.geometry = None
        self.basis = None
        self.molecular_point_group = None
        self.n_excited_states = None
        self.gs_energy_au = None
        self.excited_states = None
        self.orbs = None
        self.es_type = es_type
        with open(self.filepath) as fp:
            file = OpenFile(fp)
            self.input = extract_input_string(file)

    def read_file(self):
        with open(self.filepath) as fp:
            file = OpenFile(fp)

            read_start_time(file, self)
            self.basis = read_property(file, 'SETTING BASIS', 3)
            self.molecular_point_group = read_property(file, 'Point group', 2)
            read_structure(file, self)

            file.read_to_line_containing('PROGRAM * MULTI')

            line = file.read_to_line_containing('State symmetry')
            sym_dict = {}
            state_symmetries_read = False
            while not state_symmetries_read:
                cells = line.split()
                ix = int(cells[2])
                sym_dict[ix] = {}
                cells = file.read_to_line_containing('electrons:').split()
                sym_dict[ix]['nelec'] = int(cells[cells.index('electrons:') + 1])
                sym_dict[ix]['spin'] = cells[cells.index('Spin') + 1][9:]
                sym_dict[ix]['space sym'] = cells[cells.index('Space') + 1][9:]
                cells = file.read_to_line_containing('Number of states:').split()
                sym_dict[ix]['nstates'] = int(cells[3])
                file.read_line()
                file.read_line()
                line = file.read_line()
                if 'State symmetry' not in line:
                    state_symmetries_read = True
            space_syms = [sym_dict[i]['space sym'] for i in sym_dict.keys()]
            n_space_syms = len(list(dict.fromkeys(space_syms)))

            self.n_excited_states = 0
            self.excited_states = []
            for ix in sym_dict.keys():
                for j in range(sym_dict[ix]['nstates']):
                    cells = file.read_to_line_containing('Results for state').split()
                    state = self.es_type()
                    self.n_excited_states += 1
                    state.symmetry = cells[3][(cells[3].index('.') + 1):]
                    assert state.symmetry == sym_dict[ix]['space sym'], 'Problem reading state symmetry'
                    state.spin = sym_dict[ix]['spin']

                    cells2 = file.read_to_line_containing(f'!MCSCF STATE {cells[3]}').split()
                    if cells2[2] == '1.1':
                        self.gs_energy_au = {'CAS': float(cells2[cells2.index('Energy') + 1])}
                        state.excited_state_energy_au = float(cells2[cells2.index('Energy') + 1])
                    else:
                        state.excited_state_energy_au = float(cells2[cells2.index('Energy') + 1])

                    state.config = []
                    state.coeffs = []

                    self.excited_states.append(state)

            state_ix = 0
            for ix in sym_dict.keys():
                file.read_to_line_containing('CI')
                file.read_to_line_containing('==========================')
                file.read_line()
                ci_vec_read = False
                while not ci_vec_read:
                    line = file.read_line()
                    if line == '\n' or line == ' \n':
                        ci_vec_read = True
                    else:
                        cells = line.split()
                        for j in range(sym_dict[ix]['nstates']):
                            self.excited_states[state_ix + j].config.append(cells[0] + ' ' + cells[1])
                            self.excited_states[state_ix + j].coeffs.append(float(cells[n_space_syms + j]))
                state_ix += sym_dict[ix]['nstates']

            # if self.extension == 'caspt2':
            #    read_rs2_output(file, self, sym_dict)

            cells = file.read_to_line_containing('PROGRAMS   *').split()
            cas_ix = [i for i in range(len(cells)) if cells[i] == 'CASSCF'][0] + 1
            # if self.extension == 'caspt2':
            #    rs_ix = [i for i in range(len(cells)) if 'RS2' in cells[i]][0] + 1
            cells = file.read_line().split()
            self.cas_cpu_time = float(cells[cas_ix])
            # if self.extension == 'caspt2':
            #   self.rs2_cpu_time = float(cells[rs_ix])
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


class CasPt2OutputFile(CasOutputFile):
    def __init__(self, filepath):
        super().__init__(filepath)
        self.es_type = CasPt2ExcitedState

    def read_cas_output(self):
        self.read_file()
        return

    def read_caspt2_output(self):
        self.read_cas_output()

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            # for n_rs2_commands:
            #    read_rs2_output(file, self)


class SsSrCasPt2OutputFile(CasPt2OutputFile):
    def __init__(self, filepath):
        super().__init__(filepath)
        self.sspt2time = None
        self.mix_sets = None
        self.mixing_hams = None
        self.es_type = SsSrCasPt2ExcitedState

    def read_sssrcaspt2_output(self):
        self.read_cas_output()

        self.mix_sets = read_ms_rs2_input(self.input)
        self.mixing_hams = []

        with open(self.filepath) as fp:
            file = OpenFile(fp)

            for s in self.mix_sets:
                states = [i for i in range(len(self.excited_states))
                          if self.excited_states[i].symmetry == str(self.mix_sets[s]['symmetry'])]
                states = [states[i] for i in range(len(states))
                          if self.excited_states[states[i]].spin == self.mix_sets[s]['spin']
                          or self.excited_states[states[i]].spin == spin_dict[self.mix_sets[s]['spin']]]
                states = [j for i, j in enumerate(states)
                          if str(i + 1) + '.' + str(self.mix_sets[s]['symmetry']) in self.mix_sets[s]['states']]
                nstates = len(self.mix_sets[s]['states'])
                assert nstates == len(states)
                for j in range(nstates):
                    file.read_to_line_containing('MS-CASPT2 results after mixing')

                file.read_to_line_containing('HEFF')
                file.read_line()
                mix_ham = np.zeros((nstates, nstates))
                for i in range(nstates):
                    cells = file.read_line().split()
                    for j in range(nstates):
                        mix_ham[i, j] = float(cells[j + 1])

                file.read_to_line_containing('State        Energy      Eigenvector')
                for i in range(nstates):
                    cells = file.read_line().split()
                    self.excited_states[states[i]].sspt2_excited_state_energy_au = float(cells[1])
                    self.excited_states[states[i]].sspt2_eigenvector = np.array([
                        float(cells[2]), float(cells[3]), float(cells[4])])

        return


def read_start_time(file: OpenFile, results_object):
    line = file.read_to_line_containing('DATE')
    cells = line.split()
    results_object.start_time = {'day': cells[5][:2],
                                 'month': cells[5][3:6],
                                 'year': cells[5][-2:],
                                 'time': cells[7]}
    return


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
            results_object.geometry.append({'atom': cells[1],
                                            'xyz': np.array([float(cells[3]),
                                                             float(cells[4]),
                                                             float(cells[5])])})
    return


def read_rs2_output(file: OpenFile, results_object):
    file.read_to_line_containing('1PROGRAM * RS2')

    cells = file.read_to_line_containing('Level shift=').split()
    level_shift = float(cells[2])
    cells = file.read_to_line_containing('Reference symmetry:').split()
    ref_sym = {}
    ref_sym['space'] = cells[2]
    ref_sym['spin'] = cells[3]

    section_read = 0
    while not section_read:
        line = file.read_line()
        if '********************************************' in line:
            section_read = 1
        if 'RESULTS FOR STATE' in line:
            cells = line.split()
            sym = int(cells[3][cells[3].index('.'):])
            num = int(cells[3][:cells[3].index('.')])
            assert sym == ref_sym['space']
            states = [i for i in range(len(results_object.excited_states))
                      if results_object.excited_states[i].symmetry == sym]
            states = [states[i] for i in range(len(states))
                      if results_object.excited_states[states[i]].spin == ref_sym['spin']]
            state = states[num - 1]

            cells = file.read_to_line_containing('!RSPT2 STATE').split()
            if cells[2] == '1.1' and ref_sym['spin'] == 'Singlet':
                results_object.gs_energy_au['CASPT2'] = float(cells[4])

            results_object.excited_states[state].level_shift = level_shift
            results_object.excited_states[state].pt2_excited_state_energy_au = float(cells[4])

    return


def read_ms_rs2_input(input):
    """

    :return:

    >>> inp= '{rs2c,shift=0.3,mix=3,INIT;' + \
              'wf, charge=0,symmetry=1,spin=0;' + \
              'state,1,1}' + \
              '{rs2c,shift=0.3,mix=3;' + \
              'wf, charge=0,symmetry=1,spin=0;' + \
              'state,1,2}' + \
              '{rs2c,shift=0.3,mix=3;' + \
              'wf, charge=0,symmetry=1,spin=0;' + \
              'state,1,3}' + \
              '{rs2c,shift=0.3,mix=3,INIT;' + \
              'wf, charge=0,symmetry=1,spin=2;' + \
              'state,1,1}' + \
              '{rs2c,shift=0.3,mix=3;' + \
              'wf, charge=0,symmetry=1,spin=2;' + \
              'state,1,2}' + \
              '{rs2c,shift=0.3,mix=3;' + \
              'wf, charge=0,symmetry=1,spin=2;' + \
              'state,1,3}'
    >>> test = read_ms_rs2_input(inp)
    >>> test
    {0: {'symmetry': 1, 'spin': 0, 'states': ['1.1', '2.1', '3.1']}, 1: {'symmetry': 1, 'spin': 2, 'states': ['1.1', '2.1', '3.1']}}
    """
    segments = {}
    seg_ix = 0
    segmentation_finished = 0
    split_ix1 = input.find('rs2')
    while not segmentation_finished:
        split_ix2 = input[split_ix1 + 3:].find('rs2')
        if split_ix2 == -1:
            segments[seg_ix] = input[split_ix1:]
            segmentation_finished = 1
        else:
            segments[seg_ix] = input[split_ix1: split_ix2 + split_ix1 + 3]
        split_ix1 += split_ix2 + 3
        seg_ix += 1

    mix_ix = 0
    mix_sets = {}
    for seg in segments.keys():
        text = segments[seg]

        if 'INIT' in text:
            mix_sets[mix_ix] = {}
            text.replace(" ", "")
            j = text.index('mix')
            k = find_next_nondigit_character(text[j:])
            nstates_in_mix = int(text[j + 4: j + k + 1])

            j = text.index('symmetry')
            k = find_next_nondigit_character(text[j:])
            symmetry = text[j + 9: j + k + 1]
            mix_sets[mix_ix]['symmetry'] = int(symmetry)

            j = text.index('spin')
            k = find_next_nondigit_character(text[j:])
            mix_sets[mix_ix]['spin'] = int(text[j + 5: j + k + 1])

            mix_sets[mix_ix]['states'] = []
            for i in range(nstates_in_mix):
                text = segments[seg + i]
                text.replace(" ", "")
                j = text.index('state') + len(symmetry) + 7
                k = find_next_nondigit_character(text[j:])
                mix_sets[mix_ix]['states'].append(text[j: j + k + 1] + '.' + symmetry)
            mix_ix += 1

    return mix_sets


def extract_input_string(file: OpenFile):
    file.read_to_line_containing('***,')
    input = ''
    input_end = False
    while not input_end:
        line = file.read_line()
        if 'Checking input...' in line:
            input_end = True
        elif not line:
            break
        else:
            input += line

    return input


def read_rs2_result(results_obj: MolProOutputFile, file: OpenFile):
    res_ix = get_len_of_current_result_list(results_obj, 'CASPT2') + 1
    results_obj.results['CASPT2'][res_ix] = {}

    cells = file.read_to_line_containing('Level shift=').split()
    level_shift = float(cells[2])
    results_obj.results['CASPT2'][res_ix]['Level shift'] = level_shift
    print('\'Level shift\' added to \'CASPT2\' ' + str(res_ix) + ' results in results dictionary')

    cells = file.read_to_line_containing('Number of reference states:').split()
    rstate_ix = int(cells[6])

    cells = file.read_to_line_containing('Reference symmetry:').split()
    ref_sym = {}
    ref_sym['space'] = cells[2]
    ref_sym['spin'] = cells[3]

    if 'CAS' in results_obj.results:
        cas_keys = [i for i in results_obj.results['CAS'].keys()]
        most_recent = max(cas_keys)
        ref_states = results_obj.results['CAS'][most_recent]['Excited states']
    else:
        exit('No CAS reference results for CASPT2')

    results_obj.results['CASPT2'][res_ix]['Excited states'] = []
    results_obj.results['CASPT2'][res_ix]['Reference CAS'] = []
    section_read = 0
    escape = 0
    while not section_read:
        line = file.read_line()
        escape += 1
        if escape > 200:
            break
        if '********************************************' in line:
            section_read = 1
        if 'RESULTS FOR STATE' in line:
            pt2_state = ExcitedState()
            pt2_state.symmetry = ref_sym['space']
            pt2_state.spin = ref_sym['spin']

            cells = line.split()
            sym = cells[3][cells[3].index('.') + 1:]
            num = int(cells[3][:cells[3].index('.')])
            assert sym == ref_sym['space']
            assert num == rstate_ix

            rstates = [i for i in range(len(ref_states))
                      if ref_states[i].symmetry == sym]
            rstates = [rstates[i] for i in range(len(rstates))
                      if ref_states[rstates[i]].spin == ref_sym['spin']]
            rstate = ref_states[rstates[num - 1]]
            assert ref_sym['space'] == rstate.symmetry, 'Problem finding CAS reference results for CASPT2'
            assert ref_sym['spin'] == rstate.spin, 'Problem finding CAS reference results for CASPT2'
            results_obj.results['CASPT2'][res_ix]['Reference CAS'].append(rstate)

            cells = file.read_to_line_containing('Reference energy').split()
            ref_en = float(cells[2])
            results_obj.results['CASPT2'][res_ix]['Reference energy'] = ref_en
            print('\'Reference energy\' added to \'CASPT2\' ' + str(res_ix) + ' results in results dictionary')

            cells = file.read_to_line_containing('!RSPT2 STATE').split()
            pt2_state.excited_state_energy_au = float(cells[4])
            pt2_state.excitation_energy_eV = pt2_state.excited_state_energy_au / hartree_per_ev

            results_obj.results['CASPT2'][res_ix]['Excited states'].append(pt2_state)

    print('\'Reference CAS\' added to \'CASPT2\' ' + str(res_ix) + ' results in results dictionary')
    print('\'Excited states\' added to \'CASPT2\' ' + str(res_ix) + ' results in results dictionary')

    return



def read_cas_result(results_obj: MolProOutputFile, file: OpenFile):
    res_ix = get_len_of_current_result_list(results_obj, 'CAS') + 1
    results_obj.results['CAS'][res_ix] = {}

    file.read_to_line_containing('State symmetry')
    sym_dict = get_sym_dict(file)

    excited_states = read_cas_energies(file, sym_dict)
    results_obj.results['CAS'][res_ix]['Excited states'] = excited_states
    print('\'Excited states\' (energies) added to \'CAS\' ' + str(res_ix) + 'results in results dictionary')

    dm = read_cas_dipole_moment_matrix(file, sym_dict)
    results_obj.results['CAS'][res_ix]['Dipole Moment matrix / au'] = dm
    print('\'Dipole Moment matrix / au\' added to \'CAS\' results in results dictionary')

    read_cas_ci_vectors(file, sym_dict, results_obj.results['CAS'][res_ix]['Excited states'])
    print('\'Excited states\' (CI vectors) added to \'CAS\' ' + str(res_ix) + 'results in results dictionary')

    return


def read_cas_energies(file: OpenFile, sym_dict):
    excited_states = []
    for ix in sym_dict.keys():
        for j in range(sym_dict[ix]['nstates']):
            cells = file.read_to_line_containing('Results for state').split()
            state = ExcitedState()
            state.symmetry = cells[3][(cells[3].index('.') + 1):]
            assert state.symmetry == sym_dict[ix]['space sym'], 'Problem reading state symmetry'
            state.spin = sym_dict[ix]['spin']

            cells2 = file.read_to_line_containing(f'!MCSCF STATE {cells[3]}').split()
            state.excited_state_energy_au = float(cells2[cells2.index('Energy') + 1])

            state.config = []
            state.coeffs = []

            excited_states.append(state)

    energies = np.array([state.excited_state_energy_au for state in excited_states])
    gs_energy = np.amin(energies)
    for state in excited_states:
        state.excitation_energy_au = state.excited_state_energy_au - gs_energy
        state.excitation_energy_eV = state.excitation_energy_au / hartree_per_ev

    return excited_states


def read_cas_dipole_moment_matrix(file: OpenFile, sym_dict):
    line = file.read_to_line_containing('!MCSCF expec')
    dm = {}
    for ix in sym_dict.keys():
        dm['sym dict ' + str(ix)] = {}

    for axis in ['X', 'Y', 'Z']:
        for ix in sym_dict.keys():
            dm['sym dict ' + str(ix)][axis] = np.zeros((sym_dict[ix]['nstates'], sym_dict[ix]['nstates']))
            for i in range(sym_dict[ix]['nstates']):
                if line == '\n' or line == ' \n':
                    print('Not all expected matrix elements present in dipole moment matrix')
                    break
                for j in range(i + 1):
                    if line == '\n' or line == ' \n':
                        print('Not all expected matrix elements present in dipole moment matrix')
                        break
                    cells = line.split()
                    mat_elem = cells[2]
                    dot1 = find_next_nondigit_character(mat_elem[1:]) + 1
                    start_of_ket = int((len(mat_elem) + 1) / 2 + 2)
                    dot2 = find_next_nondigit_character(mat_elem[start_of_ket:]) + start_of_ket
                    bra = mat_elem[1:dot1]
                    ket = mat_elem[start_of_ket:dot2]
                    dm['sym dict ' + str(ix)][axis][int(bra) - 1, int(ket) - 1] = float(cells[3])
                    line = file.read_line()
        line = file.read_line()

    return dm


def read_cas_ci_vectors(file: OpenFile, sym_dict, estates):
    space_syms = [sym_dict[i]['space sym'] for i in sym_dict.keys()]
    n_space_syms = len(list(dict.fromkeys(space_syms)))

    state_ix = 0
    for ix in sym_dict.keys():
        file.read_to_line_containing('CI vector')
        file.read_to_line_containing('==========================')
        file.read_line()
        ci_vec_read = False
        while not ci_vec_read:
            line = file.read_line()
            if line == '\n' or line == ' \n':
                ci_vec_read = True
            else:
                cells = line.split()
                for j in range(sym_dict[ix]['nstates']):
                    estates[state_ix + j].config.append(cells[0] + ' ' + cells[1])
                    estates[state_ix + j].coeffs.append(float(cells[n_space_syms + j]))
        state_ix += sym_dict[ix]['nstates']
    return


def get_sym_dict(file: OpenFile):
    """Read MolPro's definitions of each symmetry index and store in a dictionary"""
    sym_dict = {}
    state_symmetries_read = False
    while not state_symmetries_read:
        cells = file.current_line.split()
        ix = int(cells[2])
        sym_dict[ix] = {}
        cells = file.read_to_line_containing('electrons:').split()
        sym_dict[ix]['nelec'] = int(cells[cells.index('electrons:') + 1])
        sym_dict[ix]['spin'] = cells[cells.index('Spin') + 1][9:]
        sym_dict[ix]['space sym'] = cells[cells.index('Space') + 1][9:]
        cells = file.read_to_line_containing('Number of states:').split()
        sym_dict[ix]['nstates'] = int(cells[3])
        file.read_line()
        file.read_line()
        line = file.read_line()
        if 'State symmetry' not in line:
            state_symmetries_read = True
    return sym_dict


class Orbital:

    def __init__(self):
        self.index = None
        self.sym_ix = None
        self.occupation = None
        self.eigenvalue = None
