#SpatialModulation.py
#Needs to be refined, in practice f_U and f_omega curves are found experimentally

import numpy as np
from scipy.integrate import simps
from scipy.constants import pi

def decreasingh(timearray, h0, h_final):
    return np.linspace(h0, h_final, len(timearray))

def mod_position(timearray, harray, f=8e7):
    return harray * np.sin(2 * pi * f * timearray)

def f_U(w, depth_array, t_array, harray, f=8e7):
    x_t = mod_position(t_array, harray, f)
    integrand = depth_array * np.exp(-2 * x_t**2 / w**2)
    f_u = simps(integrand, t_array) / np.mean(depth_array)
    return f_u

def f_omega(w, depth_array, t_array, harray, f=8e7):
    x_t = mod_position(t_array, harray, f)
    gaussian = np.exp(-2 * x_t**2 / w**2)
    second_deriv_term = (4 * x_t**2 / w**4 - 2 / w**2) * gaussian
    integrand = depth_array * second_deriv_term
    normalization = np.mean(depth_array) * (-2 / w**2)
    f_w = simps(integrand, t_array) / normalization
    return f_w




