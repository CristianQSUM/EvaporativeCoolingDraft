#Demo.py
import numpy as np
import scipy.constants as constants
from scipy.constants import pi, speed_of_light as cLight, h, atomic_mass, Boltzmann, epsilon_0, hbar
from scipy.integrate import simpson as simps
import matplotlib.pyplot as plt

#Assuming you have SimpleSpatialModulation.py, EvaporativeFunctions.py, CrossedDipoleTrap.py
from SimpleSpatialModulation import decreasingh, mod_position, f_U, f_omega
from EvaporativeFunctions import *
from CrossedDipoleTrap import CrossedDipoleTrap

a_0 = constants.physical_constants['Bohr radius'][0]
a = 98*a_0
m = 86.909180520 * atomic_mass
taulife = 26.2348E-9
naturallinewidth = 1/(2*pi*taulife)
wavelength = 1064e-9
Delta = cLight/wavelength - cLight/780e-9 #detuning
E_r = (h**2/wavelength**2)/(2*m)
k_Boltz = Boltzmann


#Set Power Arrays, Time Scale and Beam Waists
t_array = np.linspace(0, 1, 1000)
P_initial = 5  #W
P_final = 0.5  #W
P_1_array = P_initial * np.exp(np.log(P_final/P_initial) * (t_array - t_array[0]) / (t_array[-1] - t_array[0]))
P_2_array = P_1_array.copy()
W_1 = W_2 = 50e-6 #m


#Calculated Un(spatially)modulated Depth and Frequencies, power ramp still exists
unmodtrap = CrossedDipoleTrap(
    timearray=t_array,
    power1=P_1_array,
    power2=P_2_array,
    waist1=W_1,
    waist2=W_2)

unmodbeam1depth = unmodtrap.beam1depth
unmodbeam2depth = unmodtrap.beam2depth
umodtrap_depth = unmodtrap.trapdepth
unmodomega_x = unmodtrap.omegax
unmodomega_y = unmodtrap.omegay
unmodomega_z = unmodtrap.omegaz
unmodomega_bar = unmodtrap.omegabar
umodomega_bar_dot_over_omega_bar = unmodtrap.trapfrequencymodulation


#Set Spatial Modulation parameters
AOMfreq = 8e7 #AOM frequency, Hz
h_max_1 = 4*W_1 #maximum spatial modulation of beam 1, m
h_max_2 = 4*W_2 #maximum spatial modulation of beam 2, m
h_min_1 = 0 #minimum spatial modulation
h_min_2 = 0 #minimum spatial modulation
h_array_1 = decreasingh(t_array, h_max_1, h_min_1)
h_array_2 = decreasingh(t_array, h_max_2, h_min_2)
f_U_1 = f_U(W_1, unmodbeam1depth, t_array, h_array_1, f=AOMfreq)
f_U_2 = f_U(W_2, unmodbeam2depth, t_array, h_array_2, f=AOMfreq)
f_omega_1 = f_omega(W_1, unmodbeam1depth, t_array, h_array_1, f=AOMfreq)
f_omega_2 = f_omega(W_2, unmodbeam2depth, t_array, h_array_2, f=AOMfreq)

#Calculate Spatially Modulated depths and frequencies with same power ramp
modtrap = CrossedDipoleTrap(
    timearray=t_array,
    power1=P_1_array,
    power2=P_2_array,
    waist1=W_1,
    waist2=W_2,
    fU1=f_U_1,
    fU2=f_U_2,
    fw1=f_omega_1,
    fw2=f_omega_2)

mod_trap_depth = modtrap.trapdepth
mod_omega_x = modtrap.omegax
mod_omega_y = modtrap.omegay
mod_omega_z = modtrap.omegaz
mod_omega_bar = modtrap.omegabar
mod_omega_bar_dot_over_omega_bar = modtrap.trapfrequencymodulation

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs[0, 0].plot(t_array, umodtrap_depth, label='Unmodulated Trap Depth')
axs[0, 0].plot(t_array, mod_trap_depth, label='Modulated Trap Depth')
axs[0, 0].set_title('Trap Depth')
axs[0, 0].set_xlabel('Time (s)')
axs[0, 0].set_ylabel('Trap Depth (J)')
axs[0, 0].legend()
axs[0, 1].plot(t_array, unmodomega_x, label='Unmodulated Omega_x')
axs[0, 1].plot(t_array, mod_omega_x, label='Modulated Omega_x')
axs[0, 1].set_title('Omega_x')
axs[0, 1].set_xlabel('Time (s)')
axs[0, 1].set_ylabel('Frequency (Hz)')
axs[0, 1].legend()
axs[1, 0].plot(t_array, unmodomega_y, label='Unmodulated Omega_y')
axs[1, 0].plot(t_array, mod_omega_y, label='Modulated Omega_y')
axs[1, 0].set_title('Omega_y')
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].set_ylabel('Frequency (Hz)')
axs[1, 0].legend()
axs[1, 1].plot(t_array, unmodomega_z, label='Unmodulated Omega_z')
axs[1, 1].plot(t_array, mod_omega_z, label='Modulated Omega_z')
axs[1, 1].set_title('Omega_z')
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].set_ylabel('Frequency (Hz)')
axs[1, 1].legend()
plt.suptitle('Trap Depth and Frequencies')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
plt.savefig('Trap_Depth_and_Frequencies.png')



