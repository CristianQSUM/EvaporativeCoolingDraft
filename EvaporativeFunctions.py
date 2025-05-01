#EvaporativeFunctions.py
#Evaporative Cooling imports and Functions

class CrossedDipoleTrap:
    import numpy as np
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
        self.beam1depth = single_beam_depth(power1, waist1, fU1)
        self.beam2depth = single_beam_depth(power2, waist2, fU2)
        self.trapdepth = trapdepth(beam1depth, beam2depth)
        self.omegax = omega_x(beam_frequency_squared_in_mod_direction(power1, waist1, fU1), beam_frequency_squared_in_prop_direction(power2, waist2, fU2))
        self.omegay = omega_y(beam_frequency_squared_in_vert(power1, waist1, fU1), beam_frequency_squared_in_vert(power2, waist2, fU2))
        self.omegaz = omega_z(beam_frequency_squared_in_prop_direction(power1, waist1, fU1), beam_frequency_squared_in_mod_direction(power2, waist2, fU2))
        self.omegabar = geometric_mean_freq(omegax, omegay, omegaz)
        self.trapfrequencymodulation = omega_bar_dot_over_omega_bar(omegax, omegay, omegaz)
    
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

    def beam_frequency_squared_in_prop_direction(powerarray, waist, f_U_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        depth = single_beam_depth(powerarray, waist, f_U_array, atomicpolarizability)
        z_R = Rayleigh(waist)
        omega_squared = 2*depth/(mass*z_R**2)
        return(omega_squared)

    def beam_frequency_squared_in_mod_direction(powerarray, waist, f_w_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        omega_squared = (8*atomicpolarizability*powerarray*f_w_array)/(pi*mass*waist**4)
        return(omega_squared)

    def beam_frequency_squared_in_vert(powerarray, waist, f_U_array = 1.0, atomicpolarizability = alpha_natural, mass = m):
        depth = single_beam_depth(powerarray, waist, f_U_array, atomicpolarizability)
        omega_squared = 4*depth/(mass*waist**2)
        return(omega_squared)

    def omega_x(beam1modomegasquared, beam2propomegasquared):
        return(np.sqrt(beam1modomegasquared + beam2propomegasquared))

    def omega_y(beam1vertsquared, beam2vertsquared):
        return(np.sqrt(beam1vertsquared+beam2vertsquared))

    def omega_z(beam1propomegasquared, beam2modomegasquared):
        return(np.sqrt(beam1propomegasquared + beam2modomegasquared))

    def geometric_mean_freq(omegax, omegay, omegaz):
        return(np.cbrt(omegax*omegay*omegaz))

    def omega_bar_dot_over_omega_bar(omegax, omegay, omegaz):
        omegabar = geometric_mean_freq(omegax, omegay, omegaz)
        return(np.gradient(omegabar)/omegabar)

class NTdependentfunctions:
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
    
    def eta_ev(T, trap_depth, k_Boltz = Boltzmann):
        eta = trap_depth/(k_Boltz*T)
        return(eta)

    def peak_density(N, T, geometric, depth, mass = m): 
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
        n_0 = peak_density(N, T, trap_depth, geometric_frequency_array, mass)
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
        density = peak_density(N, T, geometric, depth, mass)
        approx3bodayrate = (K_3)*density**2
        return(approx3bodyrate)

    def Gamma_sc(trapdepth, detuning = Delta, Gamma = naturallinewidth):
        SpontRate = Gamma*trapdepth/(hbar*detuning)
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
        gammasc = Gamma_sc(trapdepth)
        
        eta = eta_ev(T, trap_depth)
        Efficiency = (eta + (eta-5)/(eta-4) - 3)
        dTdt = -((gammaev/3)*Effficiency - gamma3b/3 - ModulationTerm)*T + (gammasc*RecoilEnergy)/k_Boltz
        return(dTdt)
