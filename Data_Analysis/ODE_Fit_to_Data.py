import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.optimize import minimize
from scipy.integrate import odeint
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def convert_letters_to_numbers(string):

	string = str.split(string, '+ ')

	output = []

	for i in range(len(string)):

		character = string[i].strip()
		value = ord(character) - 65
		output.append(value)

	return np.asarray(output)

def interpret_single_string(string):

	reaction_k = str.split(string, ', ')
	reaction = reaction_k[0]
	k = reaction_k[1]
	
	reaction_complete = str.split(reaction, '> ')

	reactants = reaction_complete[0]
	products = reaction_complete[1]

	reac = convert_letters_to_numbers(reactants)
	prod = convert_letters_to_numbers(products)

	k_num = np.asarray([int(re.findall(r'\d+', k)[0]) - 1])

	return [reac, prod, k_num]

def ODE_string_interpreter_complete(string):

	reactions = []

	for i in string:
		converted = interpret_single_string(i)
		reactions.append(converted)

	return np.asarray(reactions)

def correct_idx(idx, shapes):

	values = []

	for counter, i in enumerate(shapes):
		value = np.ones(i) * counter
		values.append(value)

	values = np.concatenate(np.asarray(values))
	idx_corrected = []
	
	for i in idx:
		idx_corrected.append(values[i])

	return np.asarray(idx_corrected, dtype='int')

def consumption_production(matrix, mode, reaction_list_provided = False, reaction_list = None):

	number_reactants = np.amax(np.r_[np.concatenate(matrix[:,0]), np.concatenate(matrix[:,1])])

	if reaction_list_provided == False:

		reactions = []
		for i in range(number_reactants+1):
			reactions.append([])

	else:
		reactions = reaction_list

	if mode == 'consume':
		ix = 0
		k = '-k'

	if mode == 'produce':
		ix = 1
		k = 'k'

	components = np.unique(np.concatenate(matrix[:,ix]))

	idx = []
	for _ in range(len(components)):
		idx.append([])

	length = []
	for i in matrix[:,ix]:
		length.append(len(i))
	length = np.asarray(length)

	for counter, i in enumerate(components):
		idx_temp = np.where(np.concatenate(matrix[:,ix]) == i)[0]
		idx_temp = correct_idx(idx_temp, length)
		idx[counter].append(idx_temp)

	for counter, z in enumerate(idx):
		idx[counter] = [np.unique(z)]

	reax = dict(zip(components, idx))

	for i in components:
		for j in reax[i][0]:

			temp = []
			k_part = [k, matrix[j][2][0]]
			temp.append(k_part)

			for l in matrix[j][0]:
				y_part = ['y', l]
				temp.append(y_part)

			reactions[i].append(temp)

	return reactions

def ODE_interpreter(matrix):

	reactions = consumption_production(matrix, 'consume')
	reactions = consumption_production(matrix, 'produce', 'True', reactions)

	return reactions

def ODE_matrix_fit_func(k, initial_state, t, matrix, ravel = True):

	def ODE_generator(y, t):

		par = {'k': k, '-k': (-1) * k, 'y': y}
		r = []

		for component in matrix:
			r_component = 0

			for component_reaction in component:
				r_temp = 1. 

				for term in component_reaction:
					r_temp *= par[term[0]][term[1]]

				r_component += r_temp
			r.append(r_component)

		return r

	sol = odeint(ODE_generator, initial_state, t)

	if ravel == True:
		sol = np.ravel(sol)

	return sol

def residual_ODE(P, initial, t, y, matrix):
	''' Computes residual from data to be fitted and fit '''

	y_fit = ODE_matrix_fit_func(P, initial, t, matrix)
	res = np.sum((y - y_fit)**2)

	return res / 1000.

def ODE_fitting(data, t, reaction_string, sampling):

	reaction_list_formatted = ODE_string_interpreter_complete(reaction_string)
	matrix = ODE_interpreter(reaction_list_formatted)
	k_number = reaction_list_formatted[-1][-1][0] + 1

	data_ravel = np.ravel(data)

	bounds_arr = []   # bounds for rate constants in optimization

	for i in range(k_number):
		bounds_arr.append([0.0, np.inf])
	bounds_tuple = tuple(bounds_arr)

	results = []

	for _ in range(sampling):
		guess = np.random.rand(k_number)/10000.
		p = minimize(fun=residual_ODE, x0=guess, args=(data[0], t, data_ravel, matrix), bounds=bounds_tuple, method = 'L-BFGS-B')
		results.append([p.fun, p.x])

	results.sort(key=lambda x: x[0])

	p_ode = results[0][1]

	print("Residual Error:", results[0][0])
	print("Rate Constants:", p_ode)

	return p_ode, matrix

def rate_bimolecular(y, t):
	''' Defines ODEs to calculate concentration time profiles for fitting, Reaction Network:
	A + B -> C
	C -> A + B
	B + C -> D '''

	R1 = -k1 * y[0] * y[1] + k2 * y[2] * y[2]
	R2 = k3 * y[1] * y[2]

	ra = R1
	rb = R1 - R2
	rc = (-1) * R1 - R2
	rd = R2

	return [ra, rb, rc, rd]

class Kinetic_Model:

	def __init__(self, name):
		self.name = name
		self.import_xlsx()

	def import_xlsx(self):

		df = pd.read_excel('../Experimental_Data/NMR_Data/%s' % self.name, sheet_name = 0, header = 0)
		data_raw = df.to_numpy()
		data_raw = data_raw[:-1] # removing last entry due to irradiation/measurement delay

		self.t = data_raw[:,0]
		self.data = data_raw[:,1:]

	def select_normalize_data(self, data, selection):

		selection = np.asarray(selection)
		data_select = data[:,selection]
		data_norm = data_select/data_select.sum(axis = 1, keepdims = True)

		return data_norm

	def ODE_fit_to_data(self, reaction_string, sampling, selection = (2,7,10), t_full_nd = 10000):

		self.data_norm = self.select_normalize_data(self.data, selection)
		self.t_full = np.linspace(self.t[0], self.t[-1], t_full_nd)

		self.p, self.matrix = ODE_fitting(self.data_norm, self.t, reaction_string, sampling)
		self.y_solved = ODE_matrix_fit_func(self.p, self.data_norm[0], self.t_full, self.matrix, ravel = False)
		self.initial_rate = np.amax(np.gradient(self.y_solved[:,0], self.t_full))/3600

		print('Initial Rate (s^-1):', self.initial_rate)

	def plot_func(self, data, time, plot_type, show_label = False, labels = None, ax = plt):

		colors = ['blue', 'orange', 'green', 'red', 'black', 'grey']

		if labels is None:
			labels = ['A', 'B', 'C', 'D', 'E', 'F']

		if show_label == True:
			for i in range(len(data[0,:])):
				ax.plot(time, data[:,i], '%s' % plot_type, color = colors[i], markersize = 10, label = labels[i])
		else:
			for i in range(len(data[0,:])):
				ax.plot(time, data[:,i], '%s' % plot_type, color = colors[i], markersize = 10)

	def plot_results(self, labels, ax = plt):

		self.plot_func(self.data_norm, self.t, '.', show_label = True, labels = labels, ax = ax)
		self.plot_func(self.y_solved, self.t_full, '', ax = ax)

		if ax != plt:
			ax.set_xlabel('Time / h')
			ax.set_ylabel('Normalized Concentration / a.u.')
			ax.legend()
			ax.set_ylim(-0.05,1)
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	def ODE_plot_generic(self, p, reaction_string, initial_state, t):

		reaction_list_formatted = ODE_string_interpreter_complete(reaction_string)
		matrix = ODE_interpreter(reaction_list_formatted)

		y = ODE_matrix_fit_func(p, initial_state, t, matrix, ravel = False)

		self.plot_func(y, t, '')

def main():

	reaction_string = ['C + C > B, k1',
					   'C + C > A, k2',
					   'B > C, k3']

	labels = [r'Complex $ \bf{2-cis}$', r'$\bf{Oxo}$  $\bf{Dimer}$' ,r'Complex $ \bf{1}$' ]

	kinetic_model = Kinetic_Model('JS_594_Kinetic_NMR_Data.xlsx')
	kinetic_model.ODE_fit_to_data(reaction_string, 10)

	fig, ax = plt.subplots()
	kinetic_model.plot_results(labels, ax = ax)

	return kinetic_model, labels


if __name__ == '__main__':
	main()
	plt.show()



