import numpy as np 
import pandas as pd
from scipy import constants as con
from scipy.optimize import least_squares
import find_nearest as fn
import Controlled_Average as ca
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def calc_photon_energy(wavelength):
	'''Energy in eV, wavelength in nm'''

	return ((con.h*con.c)/(wavelength/1e9))/(1.602*10**(-19))

def Planck_law(x, Temp):
	'''Temp in K, wavelength (x) in nm'''

	B_lambda = ((2 * con.h * con.c**2) / (x/1.0e9)**5) * (1 / (np.exp((con.h * con.c) / ((x/1.0e9) * con.k * Temp) ) - 1))

	return B_lambda

def numerical_integration(data, points, maximum):

	x = data[:,0]
	y = data[:,1]

	idx_max = fn.find_nearest(data, maximum)[0]

	integrals = []

	for i in points:
		idx_start = fn.find_nearest(data, i)[0]
		integral = np.trapz(y[idx_start:idx_max], x[idx_start:idx_max])
		integrals.append(integral)

	integrals = np.asarray(integrals)

	return integrals

def residual_scaling_fit(p, y_a, y_b):

	res = y_a - p * y_b
	return res

def scaling_fit(a, b):
	'''Determines factor to scale b to fit a''' 

	idx = fn.find_nearest(b, a[:,0])
	b_points = b[idx]

	y_a = a[:,1]
	y_b = b_points[:,1]

	p_guess = np.ones(1)
	p = least_squares(fun=residual_scaling_fit, x0=p_guess, args=(y_a, y_b), method='lm')

	return p.x

class Absorption_Spectrum_Dual_Irradiation:

	def __init__(self, wavelengths, measured_power, lamp_temperature, ir_cutoff, rate_wavelengths):
		self.wavelengths = wavelengths
		self.measured_power = measured_power
		self.lamp_temperature = lamp_temperature
		self.ir_cutoff = ir_cutoff
		self.rate_wavelengths = rate_wavelengths
		self.wavelengths_continous = np.arange(np.amin(self.wavelengths), np.amax(self.wavelengths)+1, 1.)
		self.rate_wavelengths_continous = np.arange(np.amin(self.rate_wavelengths), np.amax(self.rate_wavelengths)+1, 1.)
		self.wavelengths_power = np.c_[self.wavelengths, self.measured_power]
		self.generate_blackbody_spectrum()

	def generate_blackbody_spectrum(self):
		'''Calculation of blackbody spectrum using provided lamp temperature and Plancks law. Also compares predicted black body
		spectrum with measured power at different wavelenghts (see check_blackbody_fit)'''

		full_wavelength_range = np.arange(100., 3000., 0.5)

		blackbody_spectrum = Planck_law(full_wavelength_range, self.lamp_temperature)
		blackbody_spectrum = np.c_[full_wavelength_range, blackbody_spectrum]

		photon_energy = calc_photon_energy(full_wavelength_range)
		blackbody_photon_spectrum = blackbody_spectrum[:,1]/photon_energy
		blackbody_photon_spectrum = np.c_[full_wavelength_range, blackbody_photon_spectrum]

		blackbody_integrated = numerical_integration(blackbody_spectrum, self.wavelengths_continous, self.ir_cutoff)
		blackbody_integrated = np.c_[self.wavelengths_continous, blackbody_integrated]

		scaling_factor = scaling_fit(self.wavelengths_power, blackbody_integrated)

		blackbody_integrated_scaled = np.c_[blackbody_integrated[:,0], blackbody_integrated[:,1] * scaling_factor]

		self.blackbody_spectrum = blackbody_spectrum
		self.blackbody_photon_spectrum = blackbody_photon_spectrum
		self.blackbody_integrated = blackbody_integrated
		self.blackbody_integrated_scaled = blackbody_integrated_scaled

	def convert_absorption_spectrum_to_rate(self, abs_spectrum):
		'''Conversion of abs_spectrum to rate for dual irradiation experiment. Calculation is done by computing the overlap between
		different intervals of the blackbody spectrum (photon count as a function of wavelength) and the provided absorption spectrum.
		The sum of the overlap is considered to be proportional to the photochemical reaction rate.
		'''

		photon_spectrum = ca.controlled_avg(self.blackbody_photon_spectrum, abs_spectrum[:,0], check = True)

		print(photon_spectrum)

		idx_max = fn.find_nearest(abs_spectrum, self.ir_cutoff)[0]

		rates = []

		for i in self.rate_wavelengths_continous:
			idx_start = fn.find_nearest(abs_spectrum, i)[0]

			absorbance_interval = abs_spectrum[:,1][idx_start:idx_max]
			photon_spectrum_interval = photon_spectrum[:,1][idx_start:idx_max]

			overlap = absorbance_interval * photon_spectrum_interval
			rates.append(np.sum(overlap))

		wavelengths_rates = np.c_[self.rate_wavelengths_continous, np.asarray(rates)]

		self.wavelengths_rates = wavelengths_rates

	def fit_theoretical_rates_to_experimental(self, exp_rates):
		'''Fitting predicted rates for abs_spectrum to experimental data for rate-wavelength dependence by using a scaling factor'''

		scaling_factor = scaling_fit(exp_rates, self.wavelengths_rates)
		rates_scaled = scaling_factor * self.wavelengths_rates[:,1]

		wavelengths_rates_scaled = np.c_[self.wavelengths_rates[:,0], rates_scaled]

		self.wavelengths_rates_scaled = wavelengths_rates_scaled

		return wavelengths_rates_scaled

	def plot_wavelength_rates(self):

		plt.plot(self.wavelengths_rates[:,0], self.wavelengths_rates[:,1])

	def check_blackbody_fit(self, ax = plt):

		ax.plot(self.blackbody_integrated_scaled[:,0], self.blackbody_integrated_scaled[:,1], linewidth = 2, color = 'green', label = 'Black-body Radiation (scaled)')
		ax.plot(self.wavelengths, self.measured_power, '.', color = 'black', markersize = 20, label = 'Datapoints')

		if ax != plt:
			ax.set_xlabel('Longpass Cut-On Wavelength / nm')
			ax.set_ylabel('Measured Power / mW')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()
