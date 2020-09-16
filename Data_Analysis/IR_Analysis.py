import numpy as np
import matplotlib.pyplot as plt
import find_nearest as fn
from scipy.signal import find_peaks
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

class Constructed_IR_Spectrum:

	def __init__(self, name, hwhm, energy_range = (0, 4000.), step_size = 0.5):
		self.name = name
		self.hwhm = hwhm
		self.data = np.loadtxt(self.name, skiprows = 1)
		self.energy = np.arange(energy_range[0], energy_range[1], step_size)
		self.generate_IR_spectrum()

	def generate_IR_spectrum(self):

		absorbance = np.zeros(len(self.energy))

		for i in self.data:

			p = np.array([self.hwhm, i[0], i[1]])
			peak = self.gaussian_hwhm(p, self.energy)

			absorbance += peak

		self.absorbance = np.asarray(absorbance)

	def gaussian_hwhm(self, p, x):

		sigma = p[0]/(np.sqrt(2*np.log(2)))

		return p[2] * np.exp(-((x - p[1])**2/(2*sigma**2)))

	def plot_spectrum(self):

		plt.plot(self.energy, self.absorbance)

class IR_Spectrum:

	def __init__(self, name, name_theoretical = None):

		self.name = name
		self.data_exp = np.loadtxt(name, delimiter = ';')

		if name_theoretical is not None:
			self.name_theoretical = name_theoretical
			self.data_theo = np.loadtxt(name_theoretical, skiprows = 1)

	def pick_peaks(self, height, prominence):

		peaks, _ = find_peaks(self.data_exp[:,1], height = height, prominence = prominence)
		self.peaks = self.data_exp[peaks]

	def plot_experimental(self, ax = plt, label_offset = 0.1):

		ax.plot(self.data_exp[:,0], self.data_exp[:,1], color = 'darkgreen', linewidth = 0.7, label = 'Experimental')

		for i in self.peaks:
			ax.annotate(np.around(i[0], 1), xy = (i[0], i[1]), xytext = (i[0], i[1] + label_offset), 
				rotation = 90, size = 9, arrowprops = dict(color = 'grey', width = 0.1, headwidth = 0.2, shrink = 0.005))

	def plot_theoretical(self, ax = plt, offset = 0.96, scaling_factor = 600., name = r'[$\bf{A}$]S$_{0}$'):

		for i in self.data_theo:
			wavenumber = i[0] * offset
			oscillator_strength = i[1] / scaling_factor

			ax.plot((wavenumber, wavenumber), (0, oscillator_strength), color = 'darkred', linewidth = 1., label = r'Theoretical (scaled, horizontally and vertically) for ' + name) 

class EPR_Spectrum():

	def __init__(self, name):
		self.data = np.loadtxt(name, skiprows = 1)

	def plot_data(self, ax = plt):

		ax.plot(self.data[:,1], self.data[:,2], color = 'black')

		if ax != plt:
			ax.set_xlabel('Magnetic Field / G')
			ax.set_ylabel('Intensity / a.u.')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

def main():

	theoretical_spectrum = Constructed_IR_Spectrum('../Computational_Data/A_Intermediate/Full_Model/F-R2-R2-Down-GPOF-Disp_IR_Data.tsv', 30.)
	#theoretical_spectrum.plot_spectrum()

	fig, ax = plt.subplots(figsize = (9, 5))

	ir_spectrum = IR_Spectrum('../Experimental_Data/IR_Data/20082401_ATR_JS615_Marx.CSV', name_theoretical = '../Computational_Data/A_Intermediate/Full_Model/F-R2-R2-Down-GPOF-Disp_IR_Data.tsv')
	ir_spectrum.pick_peaks(0., 0.02)
	ir_spectrum.plot_experimental(ax = ax)	
	ir_spectrum.plot_theoretical(ax = ax)

	ax.set_ylim(-0.03, 0.8)
	ax.set_xlim(4050., -5.)

	ax.set_xlabel(r'Wavenumber / cm$^{-1}$')
	ax.set_ylabel('Absorbance')

	handles, labels = ax.get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	ax.legend(by_label.values(), by_label.keys(), loc = 'upper left')

	return fig

def secondary():

	fig, ax = plt.subplots()
	epr = EPR_Spectrum('../Experimental_Data/EPR_Data/JS-463')
	epr.plot_data(ax = ax)

	return fig

if __name__ == '__main__':
	main()
	#secondary()
	plt.show()
