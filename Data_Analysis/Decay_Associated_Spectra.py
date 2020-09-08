import os
import numpy as np 
import matplotlib.pyplot as plt
from pathlib import Path
import re
import find_nearest as fn
import Controlled_Average as ca
from Absorption_Spectrum_to_Dual_Irradiation import scaling_fit
from UV_Vis_Spectra import import_theoretical_spectra, import_plotting_parameters
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def convert_to_OD(data):
	od = data[:,1] * 5e-3 * 0.1
	return np.c_[data[:,0], od]

def calc_SAS(gs, DAS, population):

	return (DAS - (1. - population) * gs + gs)/population

class Decay_Associated_Spectra:

	def __init__(self):
		self.import_data()
		self.theoretical = import_theoretical_spectra(return_f_triplet = True)
		self.plot_parameters = import_plotting_parameters()

	def import_data(self):

		p = Path('../Experimental_Data/Picosecond_Spectroscopy/DAS_Ru-Dihydroxokomplex/')

		data_full = {}

		for file in os.listdir('../Experimental_Data/Picosecond_Spectroscopy/DAS_Ru-Dihydroxokomplex/'):
			if file.endswith('.amplitudes'):
				exp = str.split(file, '_')[2]
				data = np.genfromtxt(p/file, delimiter='	', skip_header=4)
				data_full[exp] = data

		self.data = data_full

	def calc_difference_spectrum(self, spectrum_a, spectrum_b, population, wavelength_range):

		a = convert_to_OD(spectrum_a)
		b = convert_to_OD(spectrum_b)

		idx_a = fn.find_nearest(a, wavelength_range)
		idx_b = fn.find_nearest(b, wavelength_range)

		a_sub = a[idx_a[0]:idx_a[1]]
		b_sub = b[idx_b[0]:idx_b[1]]

		arr = {-1: a_sub, 0: b_sub}

		if np.array_equal(a_sub[:,0], b_sub[:,0]) != True:
			length_diff = (len(a_sub[:,0]) - len(b_sub[:,0])) / 1e9
			short = np.int(np.floor(length_diff))
			long_arr = -1 - short
			avg_arr = ca.controlled_avg(arr[long_arr], arr[short][:,0])
			arr[long_arr] = avg_arr
			a_sub = arr[-1]
			b_sub = arr[0]

		rest = 1. - population

		dOD = (population * b_sub[:,1] + rest * a_sub[:,1]) - a_sub[:,1]
		dmOD = 1000. * dOD

		self.theoretical_difference_spectrum = np.c_[a_sub[:,0], dmOD]

	def format_DAS_data(self, polarization):

		data = self.data[polarization]
		return np.c_[data[:,0], -1000 * data[:,2], -1000 * data[:,4]]

	def plot_DAS(self, polarization, ax = plt):

		data = self.format_DAS_data(polarization)

		#ax.plot(data[:,0], data[:,1], label = r'$\tau$ = 6 ps', color = 'seagreen')
		#ax.plot(data[:,0], data[:,2], label = r'$\tau$ = 150 ps', color = 'red')


		ax.plot(data[:,0], data[:,1], label = r'$\tau$ = 6 ps', color = 'darkgreen')
		ax.plot(data[:,0], data[:,2], label = r'$\tau$ = 150 ps', color = 'darkred')
		
		ax.plot((np.amin(data[:,0]),np.amax(data[:,0])), (0, 0), '--', color = 'grey')

		if ax != plt:
			ax.set_xlabel('Wavelength / nm')
			ax.set_ylabel(r'm$\Delta$OD')
			ax.legend()

	def plot_theoretical_difference_spectrum(self, polarization, spectrum_a, spectrum_b, name = 'None', population = 0.1, wavelength_range = (400, 730), smoothing_step_size = 20., ax = plt):
		
		das_data = self.format_DAS_data(polarization)
		das_data = ca.controlled_avg(das_data, np.arange(410, 720, smoothing_step_size))

		self.calc_difference_spectrum(self.theoretical[spectrum_a], self.theoretical[spectrum_b], population, wavelength_range)
		scaling_factor = scaling_fit(das_data, self.theoretical_difference_spectrum)

		label = self.plot_parameters[spectrum_b]['Name']
		self.plot_DAS(polarization, ax = ax)
		ax.plot(self.theoretical_difference_spectrum[:,0], scaling_factor * self.theoretical_difference_spectrum[:,1], label = r'Predicted for %s (scaled)' % label, color = 'black')
		
		ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)	
		ax.legend()

def main():
	fig, ax = plt.subplots()

	das = Decay_Associated_Spectra()
	das.plot_theoretical_difference_spectrum('MA', 'three_singlet', 'g_triplet', ax = ax) #three_singlet, g_triplet

	return fig

def secondary():
	fig, ax = plt.subplots()

	das = Decay_Associated_Spectra()
	das.plot_theoretical_difference_spectrum('MA', 'three_singlet', 'f_triplet', ax = ax)

	return fig


if __name__ == '__main__':
	#main()
	secondary()
	plt.show()

