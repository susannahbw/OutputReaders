from FileHandlingObjects import OpenFile
from MolecularPropertyObjects import ExcitedState


class TurbomoleOutputFile:

    def __init__(self, filepath):
        self.filepath = filepath
        self.results = {}

    def read_file(self):

        readfuncs_dict = {"Energy of reference wave function is": read_cc2_ref_en,
                          "| sym | multi | state |          CC2 excitation energies       |  %t1   |  %t2   |":
                              read_ricc2_excitation_energies,
                          "| sym | multi | state |          ADC(2) excitation energies    |  %t1   |  %t2   |":
                              read_ricc2_excitation_energies,
                          "| Atomic coordinate, charge and isotop information |": read_geometry}
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


def read_cc2_ref_en(results_obj: TurbomoleOutputFile, file: OpenFile):
    cells = file.current_line.split()
    results_obj.results['Reference energy / hartree'] = float(cells[6])
    print("\'Reference energy / hartree\' added to results dictionary")
    return


def read_ricc2_excitation_energies(results_obj: TurbomoleOutputFile, file: OpenFile):
    results_obj.results['Excited states'] = []
    file.read_to_line_containing('+=============================================================')
    while True:
        line = file.read_line()
        if '+=============================================================' in line:
            break
        else:
            state = ExcitedState()
            cells = line.split()
            state.excitation_energy_au = float(cells[7])
            state.excitation_energy_eV = float(cells[9])
            state.symmetry = cells[1]
            state.spin = cells[3]
            state.from_states = []
            state.to_states = []
            state.coeffs = []
            results_obj.results['Excited states'].append(state)

    nstates = len(results_obj.results['Excited states'])
    for n in range(nstates):
        cells = file.read_to_line_containing('Energy:').split()
        assert (float(cells[1]) == results_obj.results['Excited states'][n].excitation_energy_au), \
            "Problem matching states in read_cc2_excitation_energies"
        file.read_to_line_containing('occ. orb.  index spin')
        file.read_line()
        while True:
            line = file.read_line()
            if '+=======================+===' in line:
                break
            else:
                cells = line.split()
                results_obj.results['Excited states'][n].from_states.append(cells[1] + cells[2])
                results_obj.results['Excited states'][n].to_states.append(cells[5] + cells[6])
                results_obj.results['Excited states'][n].coeffs.append(float(cells[9]))

    print('\'Excited states\' added to results dictionary')
    return


def read_geometry(results_obj: TurbomoleOutputFile, file: OpenFile):
    file.read_to_line_containing('atomic coordinates')
    ix = 0
    geom_dict = {}
    while True:
        line = file.read_line()
        if len(line.strip()) == 0:
            break
        else:
            cells = line.split()
            ix += 1
            geom_dict[ix] = {'Atom type': cells[3],
                             'x': float(cells[0]),
                             'y': float(cells[1]),
                             'z': float(cells[2])}

    results_obj.results['Geometry'] = geom_dict
    print('\'Geometry\' added to results dictionary')
    return
