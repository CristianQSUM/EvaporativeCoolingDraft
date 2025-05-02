# EvaporativeCooling

Model of all-optical evaporative cooling with painted potentials

This repository contains a collection of Python scripts for modeling evaporative cooling of Rubidium-87 in a crossed optical dipole trap, following the evaporative cooling model laid out in “Rapid cooling to quantum degeneracy in dynamically shaped atom traps” Richard Roy, Alaina Green, Ryan Bowler, and Subhadeep Gupta with a few changes based on the different geometry and atomic species of our setup. A pdf with a more comprehensive review of the theory behind this code is forthcoming, a few notes on the theory and assumptions present in this code are listed in the Theory/Assumptions section of the ReadMe. 

# Repository Structure

CrossedDipoleTrap.py: Class for computing trap depth and trap frequencies. As it is currently set up the first beam propagates in the z-direction and modulates in the x-direction while the second beam propagates in the x-direction and modulates in the z-direction. 

EvaporativeFunctions.py: Contains functions dependent on atom number, $N$, temperature, $T$, and the trap depth, frequencies and modulation to determine  $\dot{N}$, $\dot{T}$, evaporation parameter $\eta$, and the various loss and heating rates, including $\Gamma_{ev}$ (the evaporation rate proportional to $\Gamma_{el}$, the elastic collision rate), $\Gamma_{3b}$, the three-body loss rate, $\Gamma_{bg}$, the rate of other background collisions

SimpleSpatialModulation.py: Functions for computing how spatial modulation changes trap depth and frequency. This is actively in the process of being refined. 

Demo.py: A sample script using the functions demonstrating trap evolution and evaporative dynamics all in SI units. 

The tunable/optimizable parameters are the power ramp of each beam, the waists of each beam, the modulation of each beam, the time taken, and the initial atom number and temperature.

After the necessary import and constants, you are able to set up an arbitrary time scale, power array and waist for each beam. None of the values set are optimized (at the moment the power array in Demo.py is an exponential decay with $P_{initial} = 5 W, P_{final} = 0.5W$ and the waists are constant and equal, it is possible to add power dependent waists). 

The unmodulated trap depths and frequencies are then calculated and stored as arrays.

You are then able to set a scanning frequency (to represent painting this must be much higher that the trap frequencies (the trap frequencies are in the kHz range)) and then set the bounds of the minimum and maximum 1-D spatial modulation. Currently my decreasingh function is linear in time, representing a tightening of the modulation, this can be adjusted in SimpleSpatialModulation.py. The modulation results in fractional reductions in the effective trap depth and in the frequencies given by the respective f_U and f_omega arrays.

Using these fractional reduction coefficients the modulated trap depths and frequencies are then calculated and stored as arrays. Bear in mind we use the same power ramp as the unmodulated case.

The Trap Depth and Frequencies are plotted as functions of time for ease of comparison.

Then using solve_ivp from scipy the Number of Atoms and Temperature are solved for in the previously given time range (which needs to match the time range of the power/frequency arrays). 

A few plots are then given.

# Theory/Assumptions

The equations of state in the high eta regime (where eta is the trap depth divided by the thermal energy) with assumptions of sufficient ergodicity are well-approximated by

$$\dot{N} = -(\Gamma_{ev} + \Gamma_{3b} + \Gamma_{bg})N$$
$$\dot{T} = -(\frac{\Gamma_{ev}}{3} (\eta + \frac{\eta - 5}{\eta - 4} - 3 ) - \frac{\Gamma_{3b}}{3} - \frac{\dot{\bar{\omega}}}{\bar{\omega}})T + \frac{\Gamma_{sc} E_r}{3 k_B}$$

Here $\Gamma_{ev}$ is the rate of evaporation, which is proportional to the peak density, the mean velocity of the atoms and the scattering cross section, with some other factors, $\Gamma_{3b}$ is the rate of three-body recombination (which we are currently approximating as a constant times the peak density squared), $\Gamma_{sc}$ is the spontantaneous scattering rate (currently in EvaporativeFunctions.py we are using a far-detuned approximation, our laser is at 1064nm whereas the relevant transition is at 780nm), $\Gamma_{bg}$ is a measured rate of background collisions, the background lifetime must be independently measured soon, currently in the code we have is fixed at 0.1 and $\bar{\omega} = (\omega_x \omega_y \omega_z)^{1/3}$, the mean geometric frequency with $\dot{\bar{\omega}}$ being the rate of change of the geometric frequency dependent on the power and spatial modulation, currently in the code this is determined using a numpy gradient method. 

The trap depth and frequencies are found similarly to the Roy et. al paper with the assumption that you are able to add the depths directly and the frequencies in quadrature which is valid with the geometry and polarizations we have planned.

There are a few more assumptions made which I will discuss in a follow-up pdf, I welcome any questions/concerns about the code at cristian.ramirez@unb.ca





