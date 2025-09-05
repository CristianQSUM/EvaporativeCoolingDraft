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
c = cLight
alpha_ground_si = 7.94e-6*h #Hz/(V/m)^2 
omega_res_THz = 2*pi*377.1074635 #D1 line (THz)
omega_res_Hz = omega_res_THz*1e12
omega_laser = 2*pi*c/wavelength
alpha_detuned = (omega_res_Hz**2 * alpha_ground_si)/(omega_res_Hz**2 - omega_laser**2)
alpha_natural = alpha_detuned/(c*epsilon_0)

def eta_ev(T, trap_depth, k_Boltz = Boltzmann):
    eta = trap_depth/(k_Boltz*T)
    return(eta)

def phase_space_density(N, T, geometric, k_Boltz = Boltzmann, h=h):
    '''
    Phase space density in the harmonic limit as in pg. 258 of O'Hara's thesis.
    https://jet.physics.ncsu.edu/theses/pdf/OHara.pdf
    '''
    traposcfreq = geometric/(2*np.pi) #Converting to Hz from rads/s, see page 100 of O'Hara's Thesis
    insideterm = (h*traposcfreq)/(k_Boltz*T)
    psd = N*np.power(insideterm, 3)
    return(psd)

def thermal_db(T, mass = m, k_Boltz = Boltzmann, hbar = hbar):
    '''
    The thermal deBroglie wavelength.
    '''
    num = 2*np.pi*np.power(hbar, 2)
    den = mass*k_Boltz*T
    db = np.sqrt(num/den)
    return(db)
    
def peak_density(N, T, geometric): 
    '''
    Under the approx the peak density is 
    Phase Space Density divided by Thermal DeBroglie Wavelength Cubed
    As in page 258
    '''
    psd = phase_space_density(N, T, geometric)
    db = thermal_db(T)
    n0 = psd/np.power(db, 3)
    return(n0)

def scattering_cross_section(scatteringlength=a):
    return(8*pi*a**2)
        
def mean_speed(T, mass=m, k_Boltz=Boltzmann):
    '''
    As in page 258
    '''
    v_bar = np.sqrt(8*k_Boltz*T/(pi*mass))
    return(v_bar)
        
def Gamma_el(N, T, geometric_frequency):
    n_0 = peak_density(N, T, geometric_frequency)
    sigma = scattering_cross_section()
    v_bar = mean_speed(T)
    return(n_0*sigma*v_bar)
    
def Gamma_ev(N, T, trap_depth, geometric_frequency):
    '''
    Valid for large eta, (eta>4), it could be eta>8, I'm reviewing the literature
    '''
    eta = eta_ev(T, trap_depth)
    evaporationrate = Gamma_el(N, T, geometric_frequency) * (eta - 4) * np.exp(-eta)
    return(evaporationrate)
        
def Gamma_3b(N, T, geometric, dparam = 1.5, L_3 = 4.3e-41, mass = m):
    '''
    The treatment for three-body loss is from the 2018 paper: https://arxiv.org/pdf/quant-ph/0602010
    The dparam = 1.5 for a harmonic trap, for anharmonic traps this must be edited
    L_3 is the value given in most experimental papers, Roy et. al made the choice of naming it K_3.
    '''
    density = peak_density(N, T, geometric)
    threebodyrate = np.power(3, -dparam)*(L_3)*density**2
    return(threebodyrate)

def Gamma_sc(trapdepth, detuning = Delta, Gamma = naturallinewidth, polarizability = alpha_detuned):
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
    The background lifetime must be independently measured, this rate is 1/backgroundlifetime.
    There is a clever way to numerically estimate this, see pages 106-107 of O'Hara's Thesis
    '''
    return(rate)
        
def N_dot(N, T, trap_depth, geometric_frequency):
    '''
    From Roy et. al
    '''
    dNdt = -(Gamma_ev(N, T, trap_depth, geometric_frequency) + Gamma_3b(N, T, geometric_frequency) + Gamma_bg())*N
    return(dNdt)

def T_dot(N, T, ModulationTerm, trap_depth, geometric_frequency, RecoilEnergy=E_r):
    '''
    From Roy et. al
    '''
    #Modulation term is omega_bar_dot_over_omega_bar
    gammaev = Gamma_ev(N, T, trap_depth, geometric_frequency)
    gamma3b = Gamma_3b(N, T, geometric_frequency)
    gammasc = Gamma_sc(trap_depth)
    eta = eta_ev(T, trap_depth)
    Efficiency = (eta + (eta-5)/(eta-4) - 3)
    dTdt = -((gammaev/3)*Efficiency - gamma3b/3 - ModulationTerm)*T + (gammasc*RecoilEnergy)/k_Boltz
    return(dTdt)
