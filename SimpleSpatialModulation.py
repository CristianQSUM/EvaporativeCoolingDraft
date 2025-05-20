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
    f_u_array = []
    for hvalue in h: #will vectorize later
        t_array = np.linspace(0, 1/modfreq, nSamples)
        x_t = mod_position(t_array, hvalue, modfreq)
        integrand = np.exp(-2 * x_t**2/w**2)
        f_u = modfreq*simps(y=integrand, x=t_array, axis = 1)
        f_u_array.append(f_u)
    return f_u_array

def f_omega(w, h, f=8e7, nSamples=100): #w is waist, h is modulation amplitude
    from scipy.integrate import simpson as simps
    t_array = np.linspace(0, 1/f, nSamples)
    f_w_array = []
    for hvalue in h:
        x_t = mod_position(t_array, hvalue, f)
        integrand = ((16 * x_t**2 - 4 * w**2)/w**4) * np.exp(-2 * x_t**2/w**2)
        average_curvature = f * simps(y=integrand, x=t_array, axis = 1)
        curvature_static = (-4/w**2)
        f_w = average_curvature/curvature_static
        f_w_array.append(f_w)
    return f_w_array




