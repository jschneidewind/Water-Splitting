import numpy as np
import matplotlib.pyplot as plt
import Liquid_Phase_O2_Analysis as lp
from Reaction_ODE_Fitting import ODE_matrix_fit_func, reaction_string_to_matrix
from utility_functions import scientific_notation, plot_func
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import odeint
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from timeit import default_timer as timer

def ODE_explicit_rate_law(p, initial_state, t, matrix, flux, cross_section, ravel = False):

	def ODE_system(y, t):
		'''ODE System assuming following reaction sequence:
		B -> A (life time p[0])
		A -- hv --> B (excitation probability p[1])
		B -- hv --> C (excitation probability p[2])
		C -- hv --> D (excitation probability p[3])
		D -- hv --> E (excitation probability p[4])
		C -> A (life time p[5])
		'''

		R1 = p[0] * y[1]
		R2 = p[1] * cross_section * flux * y[0]  
		R3 = p[2] * cross_section * flux * y[1]
		R4 = p[3] * cross_section * flux * y[2]
		R5 = p[4] * cross_section * flux * y[3]
		R6 = p[5] * y[2]

		ra = R1 - R2 + R6
		rb = -R1 + R2 - R3
		rc = R3 - R4 - R6
		rd = R4 - R5
		re = R5

		return [ra, rb, rc, rd, re]

	sol = odeint(ODE_system, initial_state, t)

	if ravel is True:
		sol = np.ravel(sol)

	return sol

def normalize_data(data):
	'''Normalizing for 550 ul of 5*10^-3 mol/l solution'''

	return data/(5e-3 * 1e6)

def absorbance_to_cross_section(A):

	return A * (1/2.61e20)

class ODE_Model:

	def __init__(self, name):
		self.name = name
		self.exps = lp.convert_xlsx_to_experiments('../Experimental_Data/Liquid_Phase_O2_Data/%s' % self.name)
		self.power_levels, self.average_data, self.time = lp.average_data(self.exps['intensity'])

		actinometry = lp.Chemical_Actinometry('../Experimental_Data/Chemical_Actinometry.xlsx', '../Experimental_Data/20120613_Lumatec2_Spektren.txt')
		self.flux_factor = actinometry.p_scaled
		self.flux_levels = self.power_levels * self.flux_factor

	def ODE_fit_function(self, p, flux_levels, cross_section, initial_state, t, ravel = True, matrix = None, idx = 2, ODE_function = 'flexible'):
		'''Solving ODE system for multiple flux levels and providing correct environment for fitting to experimental data'''

		results = []

		for flux in flux_levels:

			if ODE_function == 'flexible':
				res = ODE_matrix_fit_func(p, initial_state, t, matrix, flux = np.array([flux]), cross_section = cross_section, ravel = False)
			
			else:
				res = ODE_explicit_rate_law(p, initial_state, t, matrix, flux, cross_section[0])

			results.append(res[:,idx])

		results = np.asarray(results)

		if ravel == True:
			results = np.ravel(results)

		return results

	def residual_ODE(self, p, data, flux_levels, cross_section, initial_state, t, matrix = None, idx = 2, ODE_function = 'flexible'):
		'''Residual for fitting using ODE_fit_function'''

		y_fit = self.ODE_fit_function(p, flux_levels, cross_section, initial_state, t, matrix = matrix, idx = idx, ODE_function = ODE_function)
		res = np.sum((data - y_fit)**2)

		return res / 1000.

	def fit_ODE(self, absorbance, method, reaction_string = None, guess = np.array([1/1e-5, 1., 0.001, 0.01]), idx = 2, ODE_function = 'flexible', bounds = [[1., 1000.], [0.1, 10.], [0., 1.], [0., 1.], [0., 1.], [1., 1000.]]):
		'''Fitting experimental data using ODE_fit_function.
		Flux levels are determined by proportionality of power levels and photon flux determined in Chemical_Actinometry module
		Cross sections are calculated from Absorbance in M^-1 cm^-1
		Experimental data is normalized to fraction of initial complex amount
		
		For optimization, either minimize() or differential_evolution() is used (controlled by 'method' keyword).
		If differential_evolution is used, bounds have to be provided.
	
		ODE_function determines if one-the-fly ODE system generation is used ('flexible') or if the pre-defined ODE system
		in ODE_explicit_rate_law() is used (keyword other than 'flexible' for ODE_function).

		If 'flexible' is used, a reaction_string defining the reaction network has to be provided.

		Computationally, the explicit rate law is ca. 3-4 times faster than the on-the-fly generation. This difference becomes
		significant when a computationally expensive optimization algorithm (e.g. differential evolution) is used.

		Absorbance has to be provided as a np.array

		idx controls which column of the ODE output is used for fitting (e.g. idx = 2 for reaction component C, idx = 3 for D etc.)
		'''

		data = normalize_data(self.average_data)
		data_ravel = np.ravel(data)

		matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)

		self.matrix = matrix
		self.idx = idx

		self.cross_section = absorbance_to_cross_section(absorbance)

		self.initial_state = np.zeros(reactant_number)
		self.initial_state[0] = 1.

		bound_methods = ('L-BFGS-B', 'TNC', 'SLSQP', 'Powell', 'trust-constr')

		if method in bound_methods:
			bounds_arr = []  
			for i in range(len(guess)):
				bounds_arr.append([0.0, np.inf])
			bounds_tuple = tuple(bounds_arr)

			p = minimize(fun=self.residual_ODE, x0=guess, args=(data_ravel, self.flux_levels, self.cross_section, self.initial_state, self.time, matrix, idx, ODE_function), bounds= bounds_tuple, method = method)

		elif method == 'differential_evolution':
			p = differential_evolution(func=self.residual_ODE, bounds = bounds, args = (data_ravel, self.flux_levels, self.cross_section, self.initial_state, self.time, matrix, idx, ODE_function))
	
		else:
			p = minimize(fun=self.residual_ODE, x0=guess, args=(data_ravel, self.flux_levels, self.cross_section, self.initial_state, self.time, matrix, idx, ODE_function), method = method)

		y_solved = self.ODE_fit_function(p.x, self.flux_levels, self.cross_section, self.initial_state, self.time, ravel = False, matrix = matrix, idx = idx)

		print('residual:', p.fun)
		print('p:', p.x)
		print(p)

		self.data_normalized = data
		self.p = p.x
		self.y_solved = y_solved

	def generate_model(self, p, reaction_string, absorbance, idx = 2, ODE_function = 'flexible'):

		matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)
		
		initial_state = np.zeros(reactant_number)
		initial_state[0] = 1.

		cross_section = absorbance_to_cross_section(absorbance)

		self.y_solved = self.ODE_fit_function(p, self.flux_levels, cross_section, initial_state, self.time, ravel = False, matrix = matrix, idx = idx, ODE_function = ODE_function)
		self.data_normalized = normalize_data(self.average_data)

	def calculate_initial_rates(self, data = None, flux_levels = None):

		if data is None:
			data = self.y_solved
			flux_levels = self.flux_levels

		flux_levels = flux_levels/1e17

		rates = []

		for i in data:
			rate = np.amax(np.gradient(i))
			rates.append(rate)

		rates = np.asarray(rates)/rate

		square_x, square_y, p = lp.fit_generic(flux_levels, rates, lp.power_model)

		self.initial_rates = np.c_[flux_levels * 1e17, rates]
		self.initial_rate_power_fit = np.c_[square_x * 1e17, square_y]
		self.fitted_power = p[1] 

		return p

	def initial_rate_tau_dependence(self, absorbance, tau_values, matrix = None, idx = None):

		if matrix == None:
			matrix = self.matrix
		if idx == None:
			idx = self.idx

		flux_levels = np.linspace(np.amin(self.flux_levels), np.amax(self.flux_levels), 5)
		cross_section = absorbance_to_cross_section(absorbance)
		p = np.copy(self.p)

		fits = []

		for i in tau_values:

			lamb = 1./i
			p[0] = lamb
		
			results = self.ODE_fit_function(p, flux_levels, cross_section, self.initial_state, self.time, ravel = False, matrix = matrix, idx = idx)
			initial_rate_fit = self.calculate_initial_rates(data = results, flux_levels = flux_levels)
			fits.append(initial_rate_fit[1])

		fits = np.asarray(fits)

		self.tau_fit = np.c_[tau_values, fits]
		
	def plot_ODE_fit(self, ax = plt, show_data = True, plot_type = ''):

		labels = scientific_notation(self.power_levels[::-1] * self.flux_factor)

		if show_data is True:
			plot_func(self.data_normalized[::-1], self.time, '.', ax = ax, show_labels = False)

		plot_func(self.y_solved[::-1], self.time, plot_type, labels = labels, ax = ax, show_labels = True)

		if ax != plt:
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_{2}$ Normalized Concentration')
			ax.legend(title = r'$\bf{Photon}$ $\bf{Flux}$ $\bf{/}$ $\bf{s^{-1}}$', bbox_to_anchor = (1.0, 1.02))
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0), useMathText = True)

	def plot_tau_initial_rate(self, ax = plt):

		ax.plot(self.tau_fit[:,0], self.tau_fit[:,1], 'o-', color = 'black')
		ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
		ax.set_xscale('log')

		if ax != plt:
			ax.set_xlabel(r'$\tau$ / s')
			ax.set_ylabel(r'Power $b$ ($f(x) = ax^{b}$) describing Initial Rate/Flux Relationship')

	def plot_full_ODE_solution(self, ax = plt, flux_level = None, matrix = None, p = None, reaction_string = None, absorbance = None, initial_state = None):

		if flux_level is None:
			flux_level = self.flux_levels[-1]
		if matrix is None and reaction_string is None:
			matrix = self.matrix
		if p is None:
			p = self.p
		if initial_state is None and reaction_string is None:
			initial_state = self.initial_state

		if absorbance is None:
			cross_section = self.cross_section
		else:
			cross_section = absorbance_to_cross_section(absorbance)

		if reaction_string is not None:
			matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)
			initial_state = np.zeros(reactant_number)
			initial_state[0] = 1.

		sol = ODE_matrix_fit_func(p, initial_state, self.time, matrix, flux = np.array([flux_level]), cross_section = cross_section, ravel = False)

		plot_func(sol, self.time, '', ax = ax, transpose = True, show_labels = True)

		ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

		if ax != plt:
			ax.set_xlabel(r'Time / s')
			ax.set_ylabel(r'Normalized Concentrations')
			ax.legend()

	def plot_initial_rate_flux_dependence(self, ax = plt):

		ax.plot(self.initial_rate_power_fit[:,0], self.initial_rate_power_fit[:,1], color = 'green', label = 'Fit using $f(x) = ax^{%.1f}$' % self.fitted_power, linewidth = 2)
		ax.plot(self.initial_rates[:,0], self.initial_rates[:,1], '.', markersize = 25, color = 'black', label = 'Maximum rate')

		if ax != plt:
			ax.set_xlabel(r'Photon Flux / $ s^{-1}$')
			ax.set_ylabel(r'Maximum Rate of $O_2$ Formation / Normalized')
			ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0), useMathText = True)
			ax.ticklabel_format(axis = 'x', style = 'sci', scilimits = (0,0), useMathText = True)
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()

def main():
	''' Fitting of two-photon model to kinetic data, calculation of initial rate/flux relationship, modelling of
	initial/rate flux relationship dependence on tau. O2 is reaction component C'''

	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B > C, k3, hv1, sigma1',
				       'C > D, k4, hv1, sigma1']

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')
	
	ODE_model.fit_ODE(np.array([3000.]), 'Nelder-Mead', reaction_string = reaction_string, guess = np.array([1/1e-5, 1., 0.001, 0.01]), idx = 2)

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 
	
	ODE_model.plot_ODE_fit(ax = ax[0])

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])

	fig_init, ax_init = plt.subplots()

	ODE_model.initial_rate_tau_dependence(np.array([3000.]), np.logspace(-4, 2, 20))
	ODE_model.plot_tau_initial_rate(ax = ax_init)

	#ODE_model.plot_full_ODE_solution(ax = ax_init, flux_level = 4.15e17)

	return ODE_model, fig, ax

def Bimolecular():
	''' Fitting of bimolecular reaction model to kinetic data, calculation of initial rate/flux relationship.
		 O2 is reaction component C'''

	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B + B > C, k3',
				       'C > D, k4, hv1, sigma1']

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')
	
	p = np.array([1.23768428e+03, 2.66907120e+00, 1.39506780e-03, 2.73233123e-03])
	#ODE_model.fit_ODE(np.array([3000.]), 'Nelder-Mead', reaction_string = reaction_string, guess = np.array([1/1e-5, 1., 0.001, 0.01]), idx = 2)
	ODE_model.generate_model(p, reaction_string, np.array([3000.]), idx = 2)

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 

	ODE_model.plot_ODE_fit(ax = ax[0])

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])

	return fig, ax

def STH_estimate_B_conc():
	'''Modelling of two-photon process with life time and rate constants taken from fitting of two-photon model to data,
	but changing lifetime of B to 50 ms. Plotting of full ODE solution for Flux of 2.3e16 (solar flux from 455 - 517 nm)'''

	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B > C, k3, hv1, sigma1',
				       'C > D, k4, hv1, sigma1']

	p = np.array([20., 6.61258278e-01, 1.45155452e-03, 2.78643380e-03]) # Lifetime of B: 50 ms

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')

	fig, ax = plt.subplots(figsize = (7,5))
	fig.subplots_adjust(left = 0.2, right = 0.8)

	ODE_model.plot_full_ODE_solution(ax = ax, flux_level = 2.3e16, p = p, reaction_string = reaction_string, absorbance = np.array([3000.])) # solar flux level 2.3e16

def Three_photon_mechanism_irreversible(differential_evolution_result = False):
	'''Modelling three photon reaction, O2 is reaction component D. The second photochemical reaction intermediate C
	is formed irreversibly.
	'''
	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B > C, k3, hv1, sigma1',
				       'C > D, k4, hv1, sigma1',
				       'D > E, k5, hv1, sigma1']
		
	if differential_evolution_result is False:
		p = np.array([1.20679942e+05, 1.05657358e+00, 1.19496933e-03, 7.89827829e-03, 4.69474963e-03])

	else:
		p = np.array([7.46463674e+02, 3.92012715e-01, 6.17409650e-04, 2.75038991e-03, 1.41675329e-01])

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')
	
	#ODE_model.fit_ODE(np.array([3000.]), 'differential_evolution', reaction_string = reaction_string, guess = np.array([1/1e-5, 1., 0.001, 0.01, 0.01]), idx = 3, ODE_function = 'flexible', bounds = [[1., 1000.], [0.1, 10.], [0., 1.], [0., 1.], [0., 1.]])
	ODE_model.generate_model(p, reaction_string, np.array([3000.]), idx = 3)

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 

	ODE_model.plot_ODE_fit(ax = ax[0], show_data = True)

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])

	return fig, ax

def Three_photon_mechanism_reversible():
	'''
	Modelling of three photon reaction, O2 is reaction component D. The second photochemical reaction intermediate C reverts
	back to A. Global optimization with differential_evolution, time consuming (explicit rate law should be used).

	solution (differential evolution, explicit rate law function): p = [1.01766558e+00, 9.90130153e+00, 4.27699601e-04, 1.36969863e-02,
       2.86332333e-03, 9.00810806e+02]
	'''
	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B > C, k3, hv1, sigma1',
				       'C > D, k4, hv1, sigma1',
				       'D > E, k5, hv1, sigma1',
				       'C > A, k6']
	
	p = np.array([1.01766558e+00, 9.90130153e+00, 4.27699601e-04, 1.36969863e-02, 2.86332333e-03, 9.00810806e+02])

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')

	#ODE_model.fit_ODE(np.array([3000.]), 'differential_evolution', reaction_string = reaction_string, idx = 3, ODE_function = 'explicit', bounds = [[1., 1000.], [0.1, 10.], [0., 1.], [0., 1.], [0., 1.], [1., 1000.]])
	ODE_model.generate_model(p, reaction_string, np.array([3000.]), idx = 3)

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 

	ODE_model.plot_ODE_fit(ax = ax[0], show_data = True)

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])

	return fig, ax

def Three_photon_mechanism_reversible_cubic():
	
	reaction_string = ['B > A, k1',
				       'A > B, k2, hv1, sigma1',
				       'B > C, k3, hv1, sigma1',
				       'C > D, k4, hv1, sigma1',
				       'D > E, k5, hv1, sigma1',
				       'C > A, k6']

	p = np.array([1.14852937e+02, 9.69668195e-01, 1.26034385e-03, 7.47694029e-03, 4.77365187e-03, 1.14852937e+02])

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')

	ODE_model.generate_model(p, reaction_string, np.array([3000.]), idx = 3)

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 

	ODE_model.plot_ODE_fit(ax = ax[0], show_data = False)

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])

	return fig, ax

def bimolecular_test():

	reaction_string = ['A > B, k1, hv1, sigma1',
				       'B + B > C, k2',
				       'C > D, k3, hv1, sigma1']

	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')
	
	ODE_model.fit_ODE(np.array([3000.]), 'differential_evolution', reaction_string = reaction_string, guess = np.array([0.01, 1., 0.001, 1/1e-5]), idx = 2, ODE_function = 'flexible', bounds = [[0., 1.], [0., 1.], [0., 1.]])

	fig, ax = plt.subplots(1, 2, figsize = (9.2,5.0))
	fig.subplots_adjust(wspace = 0.67, bottom = 0.11, top = 0.7, left = 0.07, right = 0.98) 

	ODE_model.plot_ODE_fit(ax = ax[0], show_data = True)

	ODE_model.calculate_initial_rates()
	ODE_model.plot_initial_rate_flux_dependence(ax = ax[1])		





if __name__ == '__main__':
	main()
	#Bimolecular()
	#STH_estimate_B_conc()
	#Three_photon_mechanism_irreversible()
	#Three_photon_mechanism_reversible()
	#Three_photon_mechanism_reversible_cubic()
	#bimolecular_test()

	plt.show()