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


class CasExcitedState:

    def __init__(self):
        self.excited_state_energy_au = None
        self.symmetry = None
        self.spin = None
        self.config = None
        self.coeffs = None

class CasPlusExcitedState:

    def __init__(self):
        self.excited_state_energy_au = None
        self.cas_excited_state_energy_au = None
        self.symmetry = None
        self.spin = None
        self.config = None
        self.coeffs = None
        self.level_shift = None


class AdcExcitedState:

    def __init__(self):
        self.excitation_energy_eV = None
        self.excited_state_energy_au = None
        self.symmetry = None
        self.spin = None
        self.from_states = None
        self.to_states = None
        self.coeffs = None
        self.trans_dip = None
        self.osc_strength = None
        self.converged = None


class VibMode:

    def __init__(self):
        self.index = None
        self.symmetry = None
        self.freq_in_inv_cm = None
        self.red_mass_amu = None
        self.force_const_mDyne_per_A = None
        self.IR_active = None
        self.IR_intensity_KM_per_mole = None
        self.Raman_active = None
        self.norm_coord_xyz = None
        self.disp_vec_r = None
