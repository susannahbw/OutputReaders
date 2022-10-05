from MolecularPropertyObjects import ExcitedState
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
