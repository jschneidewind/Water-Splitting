import os
import logging
import sys

import numpy as np

from pymcr.mcr import McrAR
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

import find_nearest as fn
from Reaction_ODE_Fitting import ODE_conc_profile, ODE_fitting
from utility_functions import plot_func

class Time_Resolved_UV_Vis:

	def __init__(self, path = '../Experimental_Data/Absorbance_Fluorescence_Data/Time_Resolved_UV_Vis'):

		self.set_up_logger()
		self.path = path
		self.wavelength, self.time, self.data = self.import_data()

	def set_up_logger(self):

		logger = logging.getLogger('pymcr')
		logger.setLevel(logging.DEBUG)

		stdout_handler = logging.StreamHandler(stream=sys.stdout) # StdOut is a "stream"; thus, StreamHandler
		stdout_format = logging.Formatter('%(message)s')  # Setting message format, removing log level and date info
		stdout_handler.setFormatter(stdout_format)
		
		logger.addHandler(stdout_handler)

	def import_data(self):
		'''Last entry (t = 1187 min) discarded as it was a control measurement after irradiating the sample overnight'''

		data_full = []
		time_full = []

		for file in os.listdir(self.path):
			if file.endswith('abs.csv'):

				time_stamp = str.split(file, '_')[-2]
				time = str.split(time_stamp, 'm')[0]
				time = np.float(time)
				time_full.append(time)

				data = np.genfromtxt('{0}/{1}'.format(self.path, file), delimiter = ';')
				data_full.append(data[:,1])

		wavelength = data[:,0]
		data_full = np.asarray(data_full)
		time_full = np.asarray(time_full)

		idx_sort = np.argsort(time_full)

		time_full = time_full[idx_sort]
		data_full = data_full[idx_sort]

		return wavelength, time_full[:-1], data_full[:-1]

	def build_concentration_guess(self, reaction_string, rate_constants, initial_state):

		self.concentration_guess = ODE_conc_profile(rate_constants, reaction_string, initial_state, self.time)

	def perform_MCR(self, start = 263., end = 600.):
		'''self.wavelength and self.data are limited to window [start, end] for MCR analysis, resulting in 
		self.wavelength_mcr and self.data_mcr'''

		idx = fn.find_nearest(self.wavelength, (start, end))

		self.wavelength_mcr = self.wavelength[idx[0]:idx[1]] 
		self.data_mcr = self.data[:,idx[0]:idx[1]]

		self.mcrar = McrAR(max_iter=200, tol_increase=2, c_regr='NNLS', st_regr='NNLS', 
					c_constraints=[ConstraintNorm()])

		self.mcrar.fit(self.data_mcr, C=self.concentration_guess, verbose=True)

	def fit_ODE_to_MCR(self, reaction_string, sampling = 5):

		self.reaction_string_fitting = reaction_string
		self.p_ode, self.matrix = ODE_fitting(self.mcrar.C_opt_, self.time, self.reaction_string_fitting, sampling)

	def plot_data(self, ax = plt):

		cm = plt.get_cmap('gnuplot')

		for counter, i in enumerate(self.data):
			ax.plot(self.wavelength, i, label = '{0} min'.format(np.int(self.time[counter])), color = cm(1.*counter/len(self.data)))

		if ax != plt:
			ax.set_xlabel('Wavelength / nm')
			ax.set_ylabel('Absorbance')
			#ax.set_title(r'Irradiation of $ \bf{3}$ in Water ($ c =  8 \times 10^{-4}$ mol/l, $l$ = 1 cm)')
			ax.legend()

	def plot_MCR_heatmap(self, ax, fig):

		x, y = np.meshgrid(self.wavelength_mcr, self.time)

		levels = MaxNLocator(nbins=20).tick_values(0, 4)
		levels_res = MaxNLocator(nbins=20).tick_values(-0.1, 0.1)

		cmap = plt.get_cmap('seismic')

		norm = BoundaryNorm(levels, ncolors=cmap.N, clip = True)
		norm_res = BoundaryNorm(levels_res, ncolors=cmap.N, clip = True)

		im_a = ax[0].pcolormesh(x, y, self.data_mcr, cmap=cmap, norm=norm)
		im_b = ax[1].pcolormesh(x, y, self.mcrar.D_opt_, cmap=cmap, norm=norm)
		im_c = ax[2].pcolormesh(x, y, self.data_mcr - self.mcrar.D_opt_, cmap=cmap, norm=norm_res)

		im_a.set_edgecolor('face')
		im_b.set_edgecolor('face')
		im_c.set_edgecolor('face')

		fig.colorbar(im_a, ax=ax[0], label = 'Absorbance')
		fig.colorbar(im_b, ax=ax[1], label = 'Absorbance')
		fig.colorbar(im_c, ax=ax[2], label = 'Residual')

		for a in ax:
			a.set_xlabel('Wavelength / nm')
			a.set_ylabel('Time / min')

	def plot_MCR_results(self, ax, labels):

		plot_func(self.mcrar.ST_opt_.T, self.wavelength_mcr, '', ax = ax[0], show_labels = True, labels = labels, transpose = True)
		ax[0].legend()

		plot_func(self.mcrar.C_opt_, self.time, 'o-', ax = ax[1], show_labels = True, labels = labels, markersize = 5, transpose = True)
		ax[1].legend()

	def plot_ODE_fit(self, ax, labels):

		time_full = np.linspace(np.amin(self.time), np.amax(self.time), 1000.)
		fitted_profile = ODE_conc_profile(self.p_ode, self.reaction_string_fitting, self.mcrar.C_opt_[0], time_full)

		plot_func(self.mcrar.C_opt_, self.time, '.', ax = ax, show_labels = True, labels = labels, markersize = 10, transpose = True)
		plot_func(fitted_profile, time_full, '', ax = ax, show_labels = False, labels = labels, transpose = True)

		ax.set_xlabel('Time / min')
		ax.set_ylabel('Normalized Concentration')
		ax.legend()

def main():

	UVVis = Time_Resolved_UV_Vis()

	reaction_string = ['C + C > B, k1',
				       'C + C > A, k2',
				       'B > C, k3']
	rate_constants = np.array([0.002, 0.00005, 0.0000000025])
	initial_state = np.array([0., 0., 1.])

	UVVis.build_concentration_guess(reaction_string, rate_constants, initial_state)
	UVVis.perform_MCR()

	UVVis.fit_ODE_to_MCR(reaction_string)

	return UVVis

def MCR_heatmap():

	UVVis = main()

	fig, ax = plt.subplots(3, figsize = (6, 7))
	fig.subplots_adjust(left = 0.12, right = 0.95, bottom = 0.07, top = 0.95, hspace = 0.15)

	UVVis.plot_MCR_heatmap(ax = ax, fig = fig)

	ax[0].text(-0.12, 1.0, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.0, 'B', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[2].text(-0.12, 1.0, 'C', transform=ax[2].transAxes, size = 20, weight='bold')

	return fig

def MCR_results():

	UVVis = main()

	fig, ax = plt.subplots(1, 2, figsize = (8, 4))
	fig.subplots_adjust(left = 0.07, right = 0.97, bottom = 0.12, top = 0.92)

	labels = [r'Complex $ \bf{2-cis}$', r'$\bf{Oxo}$  $\bf{Dimer}$' ,r'Complex $ \bf{1}$' ]

	UVVis.plot_MCR_results(ax = ax, labels = labels)

	ax[1].grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	ax[0].set_xlabel('Wavelength / nm')
	ax[0].set_ylabel('Absorbance')

	ax[1].set_xlabel('Time / min')
	ax[1].set_ylabel('Normalized Concentration')

	ax[0].text(-0.15, 1.0, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.15, 1.0, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def UVVis_data():

	UVVis = main()

	fig, ax = plt.subplots()
	UVVis.plot_data(ax = ax)

	return fig

def ODE_fit_to_MCR():

	UVVis = main()

	labels = [r'Complex $ \bf{2-cis}$', r'$\bf{Oxo}$  $\bf{Dimer}$' ,r'Complex $ \bf{1}$' ]

	fig, ax = plt.subplots()
	UVVis.plot_ODE_fit(ax = ax, labels = labels)
	ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	return fig
	
if __name__ == '__main__':
	#main()
	MCR_results()
	plt.show()



