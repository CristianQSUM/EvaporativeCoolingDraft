#EvaporativeFunctions.py
#Evaporative Cooling imports and Functions

class CrossedDipoleTrap:
    import numpy as np
    import scipy.constants as constants
    from scipy.constants import pi, speed_of_light as cLight, h, atomic_mass, Boltzmann, epsilon_0, hbar
    
    c = cLight
    alpha_si = 7.94e-6*h
    alpha_natural = alpha_si/(c*epsilon_0)

    def __init__(self, wavelength, alpha_natural, mass, k_Boltz, scattering_length):
        self.wavelength = wavelength
        self.alpha_natural = alpha_natural
        self.mass = mass
        self.k_Boltz = k_Boltz
        self.scattering_length = scattering_length
    
    def Rayleigh(waist, l=1064e-9):
        z_r = pi*waist**2/l
        return(z_r)

    def single_beam_depth(powerarray, waist, f_U_array = 1.0, atomicpolarizability=alpha_natural):
      unmodulateddepth = (2*atomicpolatizability/pi)*(powerarray/waist**2)
      epsilon_t = unmodulateddepth*f_U_array
      return(epsilon_t)

    def trap_depth(beam1depth, beam2depth):
      U_0 = beam1depth + beam2depth
      return(U_0)

    def beam_frequency_in_prop_direction(powerarray, waist, )

    def beam_frequency_in_mod_direction()

    def 


  
