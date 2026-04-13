#SpatialModulation.py
#Needs to be refined, in practice f_U and f_omega curves are found experimentally


import numpy as np
from scipy.integrate import simpson as simps
from scipy.constants import pi

def decreasingh(timearray, h0, h_final):
    return np.linspace(h0, h_final, len(timearray))

def sine_mod_position(timearray, h, f=8e7):
    return h * np.sin(2 * pi * f * timearray)

def triangle_mod_position(timearray, h, f=8e7):
    period = 1 / f
    normalized_time = np.mod(timearray, period) / period
    triangle_wave = 2 * np.abs(2 * (normalized_time - 0.5)) - 1
    return h * triangle_wave
    
def mod_position(timearray, h, f, modtype):
    """
    Can be expanded with other forms of modulation
    """
    if modtype == 'sine':
        return sine_mod_position(timearray, h, f)
    elif modtype == 'triangle':
        return triangle_mod_position(timearray, h, f)
    else:
        raise ValueError(f"Unknown modtype '{modtype}'; use 'sine' or 'triangle'.")

def f_U(w, h, modfreq=8e7, nSamples=500, modtype='sine'):
    """
    Time-averaged depth reduction factor f_U for a scanned Gaussian beam.

    w: beam waist [m]
    h: array of modulation amplitudes [m]
    modfreq: modulation frequency [Hz]
    nSamples: number of time samples per modulation period
    modtype: 'sine' or 'triangle'
    """
    f_u_array = []
    for hvalue in h:  # can be vectorized later
        t_array = np.linspace(0, 1/modfreq, nSamples)
        x_t = mod_position(t_array, hvalue, modfreq, modtype)
        integrand = np.exp(-2 * x_t**2 / w**2)
        #average over one period; factor modfreq normalizes integral over [0, 1/modfreq]
        f_u = modfreq * simps(y=integrand, x=t_array, axis=0)
        f_u_array.append(f_u)
    return f_u_array

def f_omega(w, h, f=8e7, nSamples=500, modtype='sine'):
    """
    Time-averaged curvature factor f_omega for a scanned Gaussian beam.

    w: beam waist [m]
    h: array of modulation amplitudes [m]
    f: modulation frequency [Hz]
    nSamples: number of time samples per modulation period
    modtype: 'sine' or 'triangle'
    """
    t_array = np.linspace(0, 1/f, nSamples)
    f_w_array = []
    for hvalue in h:
        x_t = mod_position(t_array, hvalue, f, modtype)
        integrand = ((16 * x_t**2 - 4 * w**2) / w**4) * np.exp(-2 * x_t**2 / w**2)
        average_curvature = f * simps(y=integrand, x=t_array, axis=0)
        curvature_static = (-4 / w**2)
        f_w = average_curvature / curvature_static
        f_w_array.append(f_w)
    return f_w_array



