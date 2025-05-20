#SpatialModulation.py
#Needs to be refined, in practice f_U and f_omega curves are found experimentally


import numpy as np
from scipy.integrate import simpson as simps
from scipy.constants import pi

def decreasingh(timearray, h0, h_final):
    return np.linspace(h0, h_final, len(timearray))

def mod_position(timearray, h, f=8e7):
    return h * np.sin(2 * pi * f * timearray)

def f_U(w, h, modfreq=8e7, nSamples = 100): #increase nSamples for better accuracy at the cost of computation time
    from scipy.integrate import simpson as simps
    t_array = np.linspace(0, 1/modfreq, nSamples)
    x_t = mod_position(t_array, h, modfreq)
    integrand = np.exp(-2 * x_t**2/w**2)
    f_u = modfreq*simps(integrand, t_array)
    return f_u

def f_omega(w, h, f=8e7, nSamples=100): #w is waist, h is modulation amplitude
    from scipy.integrate import simpson as simps
    t_array = np.linspace(0, 1/f, nSamples)
    x_t = mod_position(t_array, harray, f)
    integrand = ((16 * x_t**2 - 4 * w**2)/w**4) * np.exp(-2 * x_t**2/w**2)
    average_curvature = f * simps(integrand, t_array)
    curvature_static = (-4/w**2)
    f_w = average_curvature/curvature_static
    return f_w




