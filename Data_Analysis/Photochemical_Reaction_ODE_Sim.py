import numpy as np
import matplotlib.pyplot as plt
import Liquid_Phase_O2_Analysis as lp
from scipy.optimize import minimize
from scipy.integrate import odeint
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def normalize_data(data):
	'''Normalizing for 550 ul of 5*10^-3 mol/l solution'''

	return data/(5e-3 * 1e6)

def absorbance_to_cross_section(A):

	return A * (1/2.61e20)

def plot_func(data, t, plot_type, labels = ['A', 'B', 'C', 'D', 'E', 'F'], ax = plt, show_labels = False, transpose = False):

	colors = ['blue', 'orange', 'green', 'red', 'black', 'grey']

	if show_labels == True:

		if transpose == False:
			for counter, i in enumerate(data):
				ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = 2, label = '{:.1e}'.format(labels[counter]))
		else:
			for counter, i in enumerate(data.T):
				ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = 2, label = '{:.1e}'.format(labels[counter]))

	else:

		if transpose == False:
			for counter, i in enumerate(data):
				ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = 2)

		else:
			for counter, i in enumerate(data.T):
				ax.plot(t, i, '%s' % plot_type, color = colors[counter], markersize = 2)

class ODE_Model:

	def __init__(self, name):
		self.name = name
		self.exps = lp.convert_xlsx_to_experiments('../Experimental_Data/Liquid_Phase_O2_Data/%s' % self.name)
		self.power_levels, self.average_data, self.time = lp.average_data(self.exps['intensity'])

		actinometry = lp.Chemical_Actinometry('../Experimental_Data/Chemical_Actinometry.xlsx', '../Experimental_Data/20120613_Lumatec2_Spektren.txt')
		self.flux_factor = actinometry.p_scaled

	def ODE_biphotonic(self, p, flux, cross_section, initial_state, t):

		def ODE_system(y, t):
			'''ODE System assuming following reaction sequence:
			B -> A (life time p[0])
			A -- hv --> B (excitation probability p[1])
			B -- hv --> C (excitation probability p[2])
			C -- hv --> D (excitation probability p[3])'''

			R1 = p[0] * y[1]
			R2 = p[1] * cross_section * flux * y[0]  
			R3 = p[2] * cross_section * flux * y[1]
			R4 = p[3] * cross_section * flux * y[2]
	
			ra = R1 - R2
			rb = -R1 + R2 - R3
			rc = R3 - R4
			rd = R4

			return [ra, rb, rc, rd]

		sol = odeint(ODE_system, initial_state, t)

		return sol

	def ODE_fit_function(self, p, flux_levels, cross_section, initial_state, t, ravel = True):
		'''Solving ODE system for multiple flux levels and providing correct environment for fitting to experimental data'''

		results = []

		for flux in flux_levels:
			res = self.ODE_biphotonic(p, flux, cross_section, initial_state, t)
			results.append(res[:,2])

		results = np.asarray(results)

		if ravel == True:
			results = np.ravel(results)

		return results

	def residual_ODE(self, p, data, flux_levels, cross_section, initial_state, t):
		'''Residual for fitting using ODE_fit_function'''

		y_fit = self.ODE_fit_function(p, flux_levels, cross_section, initial_state, t)
		res = np.sum((data - y_fit)**2)

		return res / 1000.

	def fit_ODE_biphotonic(self, absorbance, method):
		'''Fitting experimental data using ODE_fit_function.
		Flux levels are determined by proportionality of power levels and photon flux determined in Chemical_Actinometry module
		Cross sections are calculated from Absorbance in M^-1 cm^-1
		Experimental data is normalized to fraction of initial complex amount
		'''
		self.flux_levels = self.power_levels * self.flux_factor
		self.cross_section = absorbance_to_cross_section(absorbance)

		print('cross:', self.cross_section)
		self.initial_state = np.array([1., 0., 0., 0.]) 

		data = normalize_data(self.average_data)
		data_ravel = np.ravel(data)

		tau = 1e-5
		lamb = 1/tau
		guess = np.array([lamb, 1., 0.001, 0.01])


		bound_methods = ('L-BFGS-B', 'TNC', 'SLSQP', 'Powell', 'trust-constr')

		if method in bound_methods:
			bounds_arr = []  
			for i in range(len(guess)):
				bounds_arr.append([0.0, np.inf])
			bounds_arr = [[1., np.inf], [0.1, 10.], [0., 1.], [0., 1.]]
			bounds_tuple = tuple(bounds_arr)

			p = minimize(fun=self.residual_ODE, x0=guess, args=(data_ravel, self.flux_levels, self.cross_section, self.initial_state, self.time), bounds=bounds_tuple, method = method)

		else:
			p = minimize(fun=self.residual_ODE, x0=guess, args=(data_ravel, self.flux_levels, self.cross_section, self.initial_state, self.time), method = method)

		# p_copy = np.copy(p.x)
		# p_copy[0] = p_copy[0] * 100
		# p_copy[2] = p_copy[2] * 100

		y_solved = self.ODE_fit_function(p.x, self.flux_levels, self.cross_section, self.initial_state, self.time, ravel = False)

		self.data_normalized = data
		self.p = p.x
		self.y_solved = y_solved

	def calculate_initial_rates(self, data = None, flux_levels = None, ax = plt, plotting = False):

		if data is None:
			data = self.y_solved
			flux_levels = self.power_levels * self.flux_factor

		flux_levels = flux_levels/1e17

		rates = []

		for i in data:

			rate = np.amax(np.gradient(i))
			rates.append(rate)

		rates = np.asarray(rates)

		square_x, square_y, p = lp.fit_generic(flux_levels, rates, lp.power_model)

		if plotting == True:
			ax.plot(flux_levels * 1e17, rates, '.')
			ax.plot(square_x * 1e17, square_y)

		return p

	def initial_rate_tau_dependence(self, absorbance, tau_values, ax = plt):

		flux_levels = np.linspace(np.amin(self.flux_levels), np.amax(self.flux_levels), 5)
		cross_section = absorbance_to_cross_section(absorbance)
		p = np.copy(self.p)

		fits = []

		for i in tau_values:

			lamb = 1./i
			p[0] = lamb
		
			results = self.ODE_fit_function(p, flux_levels, cross_section, self.initial_state, self.time, ravel = False)
			initial_rate_fit = self.calculate_initial_rates(data = results, flux_levels = flux_levels)
			fits.append(initial_rate_fit[1])

		fits = np.asarray(fits)

		self.tau_fit = np.c_[tau_values, fits]
		
	def plot_ODE_fit(self, ax = plt):

		print(self.p)

		plot_func(self.data_normalized[::-1], self.time, '.', ax = ax, show_labels = False)
		plot_func(self.y_solved[::-1], self.time, '', labels = self.power_levels[::-1] * self.flux_factor, ax = ax, show_labels = True)

		if ax != plt:
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_{2}$ Normalized Concentration')
			ax.legend(title = r'$\bf{Photon}$ $\bf{Flux}$', bbox_to_anchor = (1.0, 1.02))
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	def plot_tau_initial_rate(self, ax = plt):

		ax.plot(self.tau_fit[:,0], self.tau_fit[:,1], 'o-', color = 'black')
		ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
		ax.set_xscale('log')


		if ax != plt:
			ax.set_xlabel(r'$\tau$ / s')
			ax.set_ylabel(r'Power $b$ ($f(x) = ax^{b}$) describing Initial Rate/Flux Relationship')

	def plot_full_ODE_solution(self, ax = plt, flux_level = None):

		if flux_level is None:
			flux_level = self.flux_levels[-1]

		sol = self.ODE_biphotonic(self.p, flux_level, self.cross_section, self.initial_state, self.time)

		plot_func(sol, self.time, '', ax = ax, transpose = True)


def main():
	ODE_model = ODE_Model('Liquid_Phase_O2_Experiments_Metadata.xlsx')
	ODE_model.fit_ODE_biphotonic(3000., 'Nelder-Mead')

	fig, ax = plt.subplots(figsize = (7,5))
	fig.subplots_adjust(left = 0.2, right = 0.8)
	ODE_model.plot_ODE_fit(ax)
	
	ODE_model.calculate_initial_rates(plotting = False)


	ODE_model.initial_rate_tau_dependence(3000., np.logspace(-4, 2, 20))
	# ODE_model.plot_tau_initial_rate(ax = ax)

	# ODE_model.plot_full_ODE_solution(ax = ax, flux_level = 4.15e17)

	#fig.savefig('200809_ODE_Model.pdf', transparent = True)

	return ODE_model


if __name__ == '__main__':
	main()
	plt.show()