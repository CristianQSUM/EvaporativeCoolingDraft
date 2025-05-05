#EvaporativeFunctions.py

import numpy as np
import scipy.constants as constants
from scipy.constants import pi, speed_of_light as cLight, h, atomic_mass, Boltzmann, epsilon_0, hbar
a_0 = constants.physical_constants['Bohr radius'][0]
a = 98*a_0
m = 86.909180520 * atomic_mass
taulife = 26.2348E-9
naturallinewidth = 1/(2*pi*taulife)
wavelength = 1064e-9
Delta = cLight/wavelength - cLight/780e-9 #detuning
E_r = (h**2/wavelength**2)/(2*m)
k_Boltz = Boltzmann

    
def eta_ev(T, trap_depth, k_Boltz = Boltzmann):
    eta = trap_depth/(k_Boltz*T)
    return(eta)

def peak_density(N, T, depth, geometric, mass = m): 
    '''
    Under the approx the peak density is 
    Phase Space Density divided by Thermal DeBroglie Wavelength Cubed
    '''
    eta = eta_ev(T, depth)
    n_0 = ((mass/(2*pi*depth))**1.5)*(geometric**3)*N*eta**1.5
    return(n_0)

def scattering_cross_section(scatteringlength=a):
    return(8*pi*a**2)
        
def mean_speed(T, mass=m, k_Boltz=Boltzmann):
    v_bar = np.sqrt(8*k_Boltz*T/(pi*mass))
    return(v_bar)
        
def Gamma_el(N, T, trap_depth, geometric_frequency, mass=m):
    n_0 = peak_density(N, T, trap_depth, geometric_frequency, mass)
    sigma = scattering_cross_section()
    v_bar = mean_speed(T, mass)
    return(n_0*sigma*v_bar)
    
def Gamma_ev(N, T, trap_depth, geometric_frequency, mass=m):
    eta = eta_ev(T, trap_depth)
    evaporationrate = Gamma_el(N, T, trap_depth, geometric_frequency, mass) * (eta - 4) * np.exp(-eta)
    return(evaporationrate)
        
def Gamma_3b(N, T, depth, geometric, mass = m, K_3 = 4.3e-41):
    '''
    The true rate is (K_3*volumeintegral(density^3))/N, K_3 has units of [m^6/s]
    '''
    density = peak_density(N, T, depth, geometric, mass)
    approx3bodyrate = (K_3)*density**2
    return(approx3bodyrate)

def Gamma_sc(trapdepth, detuning = Delta, Gamma = naturallinewidth, polarizability = 7.94e-6*h):
    #At the center, the scaterring rate varies as I(r)
    omega_L = 2*pi/(1064e-9)
    omega_0 = -detuning + omega_L
    term1 = (pi*cLight**2)/(hbar * omega_0**3)
    term2 = (omega_L/omega_0)**3
    term3 = (Gamma/(-detuning) + Gamma/(detuning - 2*omega_0))**2
    I_0 = 2*trapdepth/polarizability
    SpontRate = term1*term2*term3*I_0
    return(SpontRate)

def Gamma_bg(rate=0.1):
    '''
    The background lifetime must be independently measured, this rate is 1/backgroundlifetime
    '''
    return(rate)
        
def N_dot(N, T, trap_depth, geometric_frequency, mass=m, K_3 = 4.3e-41):
    dNdt = -(Gamma_ev(N, T, trap_depth, geometric_frequency) + Gamma_3b(N, T, trap_depth, geometric_frequency, mass, K_3) + Gamma_bg())*N
    return(dNdt)

def T_dot(N, T, ModulationTerm, trap_depth, geometric_frequency, mass=m, K_3 = 4.3e-41, RecoilEnergy=E_r):
    #Modulation term is omega_bar_dot_over_omega_bar
    gammaev = Gamma_ev(N, T, trap_depth, geometric_frequency)
    gamma3b = Gamma_3b(N, T, trap_depth, geometric_frequency, mass, K_3)
    gammasc = Gamma_sc(trap_depth)
        
    eta = eta_ev(T, trap_depth)
    Efficiency = (eta + (eta-5)/(eta-4) - 3)
    dTdt = -((gammaev/3)*Efficiency - gamma3b/3 - ModulationTerm)*T + (gammasc*RecoilEnergy)/k_Boltz
    return(dTdt)
