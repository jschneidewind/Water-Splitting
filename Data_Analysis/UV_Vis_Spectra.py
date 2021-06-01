import numpy as np
import pandas as pd
import csv as csv
from scipy import constants as con
import find_nearest as fn
from utility_functions import import_lamp_spectra, insert_image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def import_uv_vis_txt(name, return_transitions = False):
	
	spectrum = np.loadtxt(name)

	if return_transitions is True:

		transitions = []

		with open(name, newline='') as file:
		    read = csv.reader(file)

		    for line in read:
		    	string = str(line).strip('[]')
		    	string = string.strip("''")
		    	split = str.split(string, '   ')
		    	
		    	if split[0] is '#':
		    		split = str.split(string, ' ')
		    		string = list(filter(None, split))

		    		if string[1] is not 'X':
		    			transitions.append([float(string[1]), float(string[2])])

		transitions = np.asarray(transitions)

		return spectrum, transitions
	
	else:
		return spectrum

def import_theoretical_spectra(return_transitions = False, return_f_triplet = False):

	if return_f_triplet is False:

		theoretical = {'three_singlet': import_uv_vis_txt('../Computational_Data/A_Intermediate/Mono/3-Full-GP-Disp-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'f_singlet': import_uv_vis_txt('../Computational_Data/A_Intermediate/Full_Model/F-R2-R2-Down-GP-Disp-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'g_r2_r2': import_uv_vis_txt('../Computational_Data/B_Intermediate/Full_Model/G-Cis-R2-R2-Down-T-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_triplet': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Full_Model/G-R2-R2-Down-GP-T-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_mono': import_uv_vis_txt('../Computational_Data/B_Intermediate/Mono/G-Cis-Mono-D-0-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_trans_mono': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Mono/G-Mono-D-0-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'three_hcis': import_uv_vis_txt('../Computational_Data/B_Intermediate/Mono/3-HCis-GP-D-0-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'three_trans_triplet': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Mono/3-HTrans-GP-D-0-TD-SMD_uvvis_full.txt', return_transitions = return_transitions),
			           '2-trans-B': import_uv_vis_txt('../Computational_Data/F_Intermediate/Trans-B/2-Trans-B-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
			           '2-trans-A': import_uv_vis_txt('../Computational_Data/F_Intermediate/Trans-A/2-Trans-A-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions)
			           }
		
		return theoretical

	else:
		theoretical = {'three_singlet': import_uv_vis_txt('../Computational_Data/A_Intermediate/Mono/3-Full-GP-Disp-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'f_singlet': import_uv_vis_txt('../Computational_Data/A_Intermediate/Full_Model/F-R2-R2-Down-GP-Disp-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'g_r2_r2': import_uv_vis_txt('../Computational_Data/B_Intermediate/Full_Model/G-Cis-R2-R2-Down-T-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_triplet': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Full_Model/G-R2-R2-Down-GP-T-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_mono': import_uv_vis_txt('../Computational_Data/B_Intermediate/Mono/G-Cis-Mono-D-0-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'g_trans_mono': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Mono/G-Mono-D-0-GP-TD-SMD_uvvis.txt', return_transitions = return_transitions),
					   'three_hcis': import_uv_vis_txt('../Computational_Data/B_Intermediate/Mono/3-HCis-GP-D-0-TD-SMD-NStates25_uvvis.txt', return_transitions = return_transitions),
					   'three_trans_triplet': import_uv_vis_txt('../Computational_Data/B_Intermediate/Trans_Isomer/Mono/3-HTrans-GP-D-0-TD-SMD_uvvis_full.txt', return_transitions = return_transitions),
					   'f_triplet': import_uv_vis_txt('../Computational_Data/A_Intermediate/Full_Model/Triplet/F-R2-R2-Down-GP-T-TD-SMD_uvvis.txt', return_transitions = return_transitions)
			           }
		
		return theoretical

def import_plotting_parameters():

	df = pd.read_excel('../Computational_Data/Misc/UV_Vis_Plots.xlsx', index_col = 0, dtype = object)
	data = df.to_dict('index')

	return data

def calc_molar_attenuation(absorbance, concentration, pathlength):

	return absorbance/(concentration*pathlength)

def import_convert_UV_Vis(name = '/UV_Vis_Complex_1/JS_594_40ms_100acc.csv', concentration = 4e-5, pathlength = 1., delimiter = ';'):

	data = np.loadtxt('../Experimental_Data/Absorbance_Fluorescence_Data/' + name, delimiter = delimiter)

	wavelength = data[:,0]
	molar_attenuation = calc_molar_attenuation(data[:,1], concentration, pathlength)

	return np.c_[wavelength, molar_attenuation]

def import_fluorescence(name):

	data = np.genfromtxt('../Experimental_Data/Absorbance_Fluorescence_Data/' + name, delimiter = ',', skip_header = 2, skip_footer = 38, usecols = (0,1))
	return data

class Experimental_UV_Vis_Fluorescence:

	def __init__(self, fluorescence_name, fluorescence_control):
		self.data = import_convert_UV_Vis()
		self.fluorescence = import_fluorescence(fluorescence_name)
		self.fluorescence_control = import_fluorescence(fluorescence_control)

	def plot_uv_vis_fluorescence(self, plotting_range = [300, 900], ax = plt, fluorescence_range = [457., 691.]):

		idx_uvvis = fn.find_nearest(self.data, plotting_range)
		self.data_plot = self.data[idx_uvvis[0]:idx_uvvis[1]]

		idx_fluorescence = fn.find_nearest(self.fluorescence, fluorescence_range)
		self.fluorescence_plot = self.fluorescence[idx_fluorescence[0]:idx_fluorescence[1]]

		plot_a = ax.plot(self.data_plot[:,0], self.data_plot[:,1], label = 'Absorbance', color = 'darkgreen', linewidth = 2.)
		
		axb = ax.twinx()
		axb.spines['right'].set_position(('axes', 1.2))

		ax.set_ylabel(r'$ \epsilon $ / $M^{-1} cm^{-1}$')
		axb.set_ylabel('Emission / a.u.')

		plot_b = axb.plot(self.fluorescence_plot[:,0], self.fluorescence_plot[:,1], label = 'Fluorescence\n(Excitation: 370 nm)', color = 'darkred', linewidth = 2.)

		return plot_a, plot_b, ax, axb

class Theoretical_UV_Vis_Spectrum:

	def __init__(self, spectra, name, plot_parameter):
		self.all_spectra = spectra
		self.name = name
		self.spectrum = spectra[name][0]
		self.transitions = spectra[name][1]
		self.plot_parameter = plot_parameter

	def plot_image(self, x, y, ax):

		image = self.plot_parameter[self.name]['Image']

		if image != 'No':

			vert_off = self.plot_parameter[self.name]['Vertical_Offset']
			horiz_off = self.plot_parameter[self.name]['Horizontal_Offset']

			img = mpimg.imread('../Computational_Data/Images/%s' % image)
			imagebox = OffsetImage(img, zoom = self.plot_parameter[self.name]['Zoom'])

			ab = AnnotationBbox(imagebox, (x + horiz_off, y + vert_off), frameon = False)
			ax.add_artist(ab)

	def plot_uv_vis(self, ax = plt, plotting_range = [300, 900], color = 'black', label = None, plotting_image = True, scaling_factor = 1., linewidth = 1.5, oscillator_ylim = False):

		if label is None:
			label = self.plot_parameter[self.name]['Name']

		idx = fn.find_nearest(self.spectrum, plotting_range)
		self.spectrum_plot = self.spectrum[idx[0]:idx[1]]
		self.spectrum_plot[:,1] = self.spectrum_plot[:,1] * scaling_factor

		plot_a = ax.plot(self.spectrum_plot[:,0], self.spectrum_plot[:,1], label = label, color = color, linewidth = linewidth)

		axb = ax.twinx()

		if oscillator_ylim is True:
			axb.set_ylim(-0.01, 0.2)

		for i in self.transitions:
			if i[0] < plotting_range[1] and i[0] > plotting_range[0]:
				axb.plot((i[0], i[0]), (0, i[1]), color = color, linewidth = linewidth)

		if ax != plt:
			ax.set_xlabel('Wavelength / nm')
			ax.set_ylabel(r'$\epsilon$ / a.u.')
			axb.set_ylabel(r'Oscillator Strength')
			ax.legend(fontsize = 'small', loc = 'upper right')

			if plotting_image == True:
				y_mid = ax.get_ylim()[1]/2
				self.plot_image(plotting_range[1], y_mid, ax = ax)

		return plot_a, axb

class Constructed_UV_Vis_Spectrum:

	def __init__(self, name, hwhm):
		self.name = name
		self.hwhm = hwhm

		df = pd.read_excel(name)
		data = df.to_numpy()

		self.transition_energies = data[:,0]
		self.oscillator_strengths = data[:,1]

		self.generate_uv_vis()

	def generate_uv_vis(self):
		'''transition energy in nm, half width half maximum in eV'''

		transition_energy = self.convert_nm_to_eV(self.transition_energies)

		energy_range_nm = np.arange(200, 2000, 0.5)
		energy_range_eV = self.convert_nm_to_eV(energy_range_nm)

		absorbance = np.zeros(len(energy_range_eV))

		for i in zip(transition_energy, self.oscillator_strengths):

			p = np.array([self.hwhm, i[0], i[1]])
			peak = self.gaussian_hwhm(p, energy_range_eV)

			absorbance += peak

		self.uv_vis = np.c_[energy_range_nm, absorbance]

	def gaussian_hwhm(self, p, x):

		sigma = p[0]/(np.sqrt(2*np.log(2)))

		return p[2] * np.exp(-((x - p[1])**2/(2*sigma**2)))

	def convert_nm_to_eV(self, wavelength):
		'''Energy in eV, wavelength in nm'''

		return ((con.h*con.c)/(wavelength/1e9))/(1.602*10**(-19))

	def convert_eV_to_nm(self, eV):

		return (con.h*con.c*1e9) / (eV * 1.602e-19)

	def plot_uv_vis(self, ax = plt):

		ax.plot(self.uv_vis[:,0], self.uv_vis[:,1])

class Light_Source_Spectrum:

	def __init__(self, name):
		self.name = name
		self.spectra_320_400, self.spectra_320_500 = import_lamp_spectra(self.name)

	def plot_spectra(self, ax = plt):

		ax.plot(self.spectra_320_400[:,0], self.spectra_320_400[:,1], label = '320 - 400 nm Filter')
		ax.plot(self.spectra_320_500[:,0], self.spectra_320_500[:,1] * 2, label = '320 - 500 nm Filter (scaled)')

		if ax != plt:
			ax.set_xlabel('Wavelength / nm')
			ax.set_ylabel(r'Power / $\mu$W cm$^{-2}$ nm$^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()

def main():

	theoretical_spectra = import_theoretical_spectra(return_transitions = True)
	plotting_parameter = import_plotting_parameters()

	vertical_plots = 4
	horizontal_plots = 2

	colors = {0: 'darkblue', 1: 'darkblue', 2: 'darkgreen', 3: 'darkred', 4: 'darkgreen', 5: 'darkred', 6: 'darkgreen', 7: 'darkred'}

	fig, ax = plt.subplots(vertical_plots, horizontal_plots, figsize = (8,10))
	fig.subplots_adjust(wspace = 0.55, hspace = 0.3, top = 0.95, bottom = 0.05)
	
	for counter, i in enumerate(theoretical_spectra):

		spectrum = Theoretical_UV_Vis_Spectrum(theoretical_spectra, i, plotting_parameter)
		pos_a = int(np.floor(counter / horizontal_plots))
		pos_b = counter - pos_a * horizontal_plots

		spectrum.plot_uv_vis(ax = ax[pos_a][pos_b], color = colors[counter])
		ax[pos_a][pos_b].text(-0.07, 1.05, '%s' % chr(65+counter), transform=ax[pos_a][pos_b].transAxes, size = 20, weight='bold')

	return fig

def secondary():

	theoretical_spectra = import_theoretical_spectra(return_transitions = True)
	plotting_parameter = import_plotting_parameters()

	fig, ax = plt.subplots()
	fig.subplots_adjust(right = 0.75, top = 0.95)

	theoretical = Theoretical_UV_Vis_Spectrum(theoretical_spectra, 'f_singlet', plotting_parameter)
	plot_a, ax_a = theoretical.plot_uv_vis(ax = ax, scaling_factor = 0.3, plotting_image = False, linewidth = 1., oscillator_ylim = True, label = 'Predicted Absorbance for\n' + r'[$\bf{A}$]S$_{0}$ (scaled)', color = 'black')

	experimental = Experimental_UV_Vis_Fluorescence('Fluorescence_Complex_1/JS591_exc370_20nm_slit_each_slow.csv', 'Fluorescence_Water/Wasser_exc370_20nm_slit_each_slow.csv')
	plot_b, plot_c, ax_b, ax_c = experimental.plot_uv_vis_fluorescence(ax = ax)

	plots = plot_a + plot_b + plot_c
	labels = [p.get_label() for p in plots]
	colors = [p.get_color() for p in plots]

	ax_b.yaxis.label.set_color(colors[1])
	ax_b.tick_params(axis = 'y', colors = colors[1])

	ax_c.yaxis.label.set_color(colors[2])
	ax_c.tick_params(axis = 'y', colors = colors[2])

	ax.legend(plots, labels)

	return fig, ax_a

def tertiary():

	fig, ax = plt.subplots()

	lumatec = Light_Source_Spectrum('../Experimental_Data/20120613_Lumatec2_Spektren.txt')
	lumatec.plot_spectra(ax = ax)

	return fig

def two_trans_uv_vis():

	theoretical_spectra = import_theoretical_spectra(return_transitions = True)
	plotting_parameter = import_plotting_parameters()


	fig, ax = plt.subplots()
	fig.subplots_adjust(right = 0.75)

	plotting_range = [230, 800]

	theoretical = Theoretical_UV_Vis_Spectrum(theoretical_spectra, '2-trans-B', plotting_parameter)
	plot_a, ax_a = theoretical.plot_uv_vis(ax = ax, scaling_factor = 0.4, plotting_image = False, 
										linewidth = 1., oscillator_ylim = True, 
										label = 'Predicted Absorbance for\n' + r'[$\bf{F-Trans-Up}$]S$_{0}$ (scaled)', 
										color = 'black', plotting_range = plotting_range)

	Milstein_2009_data = import_convert_UV_Vis(name = 'UV_Vis_Complex_2_Milstein2009/Milstein_2009_Complex_2_UV_Vis.txt',
										concentration = 4.26e-4, pathlength = 1, delimiter = ' ')

	idx_uvvis = fn.find_nearest(Milstein_2009_data, plotting_range)
	data_plot = Milstein_2009_data[idx_uvvis[0]:idx_uvvis[1]]

	plot_b = ax.plot(data_plot[:,0], data_plot[:,1], label = r'Absorbance $\bf{2-trans}$ (Milstein, 2009)', color = 'darkgreen', linewidth = 2.)


	plots = plot_a + plot_b
	labels = [p.get_label() for p in plots]
	colors = [p.get_color() for p in plots]

	ax.set_ylabel(r'$ \epsilon $ / $M^{-1} cm^{-1}$')

	ax.yaxis.label.set_color(colors[1])
	ax.tick_params(axis = 'y', colors = colors[1])

	insert_image('../Computational_Data/F_Intermediate/Trans-B/2-Full-Trans-B-GPNMR-NoSpin.png', 0.63, 0.5, 0.22, ax)

	ax.legend(plots, labels)

	return fig



if __name__ == '__main__':
	#main()
	#secondary()
	#tertiary()
	two_trans_uv_vis()
	plt.show()



