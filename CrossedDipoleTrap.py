#CrossedDipoleTrap.py
import numpy as np

class CrossedDipoleTrap:
    import scipy.constants as constants
    from scipy.constants import pi, speed_of_light as cLight, h, atomic_mass, Boltzmann, epsilon_0, hbar
    
    c = cLight
    alpha_si = 7.94e-6*h
    alpha_natural = alpha_si/(c*epsilon_0)
    m = 86.909180520 * atomic_mass
    a_0 = constants.physical_constants['Bohr radius'][0]
    a = 98*a_0


    def __init__(self, timearray, power1, power2, waist1, waist2, fU1 = 1.0, fU2 = 1.0, fw1 = 1.0, fw2 = 1.0,
                 wavelength = 1064e-9, alpha_natural = alpha_natural, mass = m, k_Boltz = Boltzmann, scattering_length = a):
        self.wavelength = wavelength
        self.alpha_natural = alpha_natural
        self.mass = mass
        self.k_Boltz = k_Boltz
        self.scattering_length = scattering_length
        self.timearray = timearray
        self.power1 = power1
        self.power2 = power2
        self.waist1 = waist1
        self.waist2 = waist2
        self.beam1depth = self.single_beam_depth(power1, waist1, fU1)
        self.beam2depth = self.single_beam_depth(power2, waist2, fU2)
        self.trapdepth = self.trap_depth(self.beam1depth, self.beam2depth)
        self.omegax = self.omega_x(self.beam_frequency_squared_in_mod_direction(power1, waist1, fU1), self.beam_frequency_squared_in_prop_direction(power2, waist2, fU2))
        self.omegay = self.omega_y(self.beam_frequency_squared_in_vert(power1, waist1, fU1), self.beam_frequency_squared_in_vert(power2, waist2, fU2))
        self.omegaz = self.omega_z(self.beam_frequency_squared_in_prop_direction(power1, waist1, fU1), self.beam_frequency_squared_in_mod_direction(power2, waist2, fU2))
        self.omegabar = self.geometric_mean_freq(self.omegax, self.omegay, self.omegaz)
        self.trapfrequencymodulation = self.omega_bar_dot_over_omega_bar(self.omegax, self.omegay, self.omegaz)
    
    def single_beam_depth(self, powerarray, waist, f_U_array = 1.0, atomicpolarizability=alpha_natural):
        from scipy.constants import pi
        unmodulateddepth = (2*atomicpolarizability/pi)*(powerarray/waist**2)
        epsilon_t = unmodulateddepth*f_U_array
        return(epsilon_t)
    
    def Rayleigh(self, waist, l=1064e-9):
        from scipy.constants import pi
        z_r = pi*waist**2/l
        return(z_r)

    def trap_depth(self, beam1depth, beam2depth):
      U_0 = beam1depth + beam2depth
      return(U_0)

    def beam_frequency_squared_in_prop_direction(self, powerarray, waist, f_U_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        depth = self.single_beam_depth(powerarray, waist, f_U_array, atomicpolarizability)
        z_R = self.Rayleigh(waist)
        omega_squared = 2*depth/(mass*z_R**2)
        return(omega_squared)

    def beam_frequency_squared_in_mod_direction(self, powerarray, waist, f_w_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        from scipy.constants import pi
        omega_squared = (8*atomicpolarizability*powerarray*f_w_array)/(pi*mass*waist**4)
        return(omega_squared)

    def beam_frequency_squared_in_vert(self, powerarray, waist, f_U_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        depth = self.single_beam_depth(powerarray, waist, f_U_array, atomicpolarizability)
        omega_squared = 4*depth/(mass*waist**2)
        return(omega_squared)

    def omega_x(self, beam1modomegasquared, beam2propomegasquared):
        return(np.sqrt(beam1modomegasquared + beam2propomegasquared))

    def omega_y(self, beam1vertsquared, beam2vertsquared):
        return(np.sqrt(beam1vertsquared+beam2vertsquared))

    def omega_z(self, beam1propomegasquared, beam2modomegasquared):
        return(np.sqrt(beam1propomegasquared + beam2modomegasquared))

    def geometric_mean_freq(self, omegax, omegay, omegaz):
        return(np.cbrt(omegax*omegay*omegaz))

    def omega_bar_dot_over_omega_bar(self, omegax, omegay, omegaz):
        omegabar = self.geometric_mean_freq(omegax, omegay, omegaz)
        return(np.gradient(omegabar)/omegabar)
        
