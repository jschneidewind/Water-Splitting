import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.optimize import minimize
from scipy.integrate import odeint
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from timeit import default_timer as timer
from utility_functions import plot_func
from numba import njit
from numba.typed import List
import pprint

def convert_letters_to_numbers(string):

	string = str.split(string, '+')

	output = []

	for item in string:
		character = item.strip(' ')
		value = ord(character) - 65
		output.append(value)

	return np.asarray(output)

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

def interpret_single_string(string):

	reaction_and_non_reactant_part = str.split(string, ',')
	reaction = reaction_and_non_reactant_part[0]
	non_reactant_part = reaction_and_non_reactant_part[1:]

	reactants, products = str.split(reaction, '>')

	non_reactant_components = []

	for item in non_reactant_part:
		item = item.strip(' ')

		component_type = re.sub(r'\d+', '', item)
		index = int(re.findall(r'\d+', item)[0]) - 1 # converting from 1 to 0 indexing

		non_reactant_components.append([component_type, index])

	return [convert_letters_to_numbers(reactants), 
			convert_letters_to_numbers(products)], non_reactant_components

def ODE_string_interpreter_complete(string):

	reactions = []
	non_reactant_components = {}

	for counter, i in enumerate(string):
		reaction, non_reactant_component = interpret_single_string(i)
		reactions.append(reaction)

		non_reactant_components[counter] = non_reactant_component

	return np.asarray(reactions), non_reactant_components

def consumption_production(matrix, non_reactant_components, mode, reaction_list_provided = False, reaction_list = None):

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
	for _ in components:
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

	mapping = {'y': 0, 'k': 1, '-k': 2, 'hv': 3, 'sigma': 4}

	for i in components:

		for j in reax[i][0]:
			temp = []

			for component_type, component_index in non_reactant_components[j]:
				if component_type == 'k':
					component_type = k

				part = [mapping[component_type], component_index] # component_types are converted to intergers based on mapping dict
				temp.append(part)

			for l in matrix[j][0]:
				y_part = [0, l]  # component type y is assigned integer 0
				temp.append(y_part)

			reactions[i].append(temp)

	return reactions

def ODE_interpreter(matrix, non_reactant_components):

	reactions = consumption_production(matrix, non_reactant_components, 'consume')
	reactions = consumption_production(matrix, non_reactant_components, 'produce', 'True', reactions)

	return reactions

def reaction_string_to_matrix(string):

	reax, non_reactant_components = ODE_string_interpreter_complete(string)

	reactant_number = len(np.unique(np.r_[np.concatenate(reax[:,0]), np.concatenate(reax[:,1])]))
	
	k_indices = []

	for key, value in non_reactant_components.items():
		for item in value:
			if 'k' in item:
				k_indices.append(item[1])

	k_number = len(np.unique(np.asarray(k_indices)))

	matrix = ODE_interpreter(reax, non_reactant_components)

	return matrix, k_number, reactant_number

def convert_nested_list_to_numba(lst):

	numba_matrix = List()

	for element in lst:
		if isinstance(element, list):
			numba_matrix.append(convert_nested_list_to_numba(element))
		else:
			numba_matrix.append(element)

	return numba_matrix

def reaction_string_to_numba_matrix(string):

	reaction_matrix, k_number, reactant_number = reaction_string_to_matrix(string)

	return convert_nested_list_to_numba(reaction_matrix), k_number, reactant_number

def ODE_generator(y, t, output, matrix, k, flux, cross_section):
	'''Generation of ODE system based on provided matrix (list or numba.typedList)
	'''
	par = {0: y, 1: k, 2: (-1) * k, 3: flux, 4: cross_section}

	for counter, component in enumerate(matrix):
		r_component = 0

		for component_reaction in component:
			r_temp = 1. 

			for term in component_reaction:
				r_temp *= par[term[0]][term[1]]

			r_component += r_temp

		output[counter] = r_component

	return output

ODE_generator_numba = njit(ODE_generator, nogil=True) # jit compilation of ODE_generator

def ODE_matrix_fit_func(k, initial_state, t, matrix, ravel = True, 
			flux = np.array([1e18]), cross_section = np.array([2e-18, 4e-18, 8e-18]), numba = False):
	'''Numerically solving ODE system generated by ODE_generator using SciPy's odeint function.
	If numba is True, jit compiled ODE_generator is requested. In this case, "matrix" has to be a
	numba.TypedList (generated by reaction_string_to_numba_matrix)
	'''

	if numba is True:
		ODE_function = ODE_generator_numba
	else:
		ODE_function = ODE_generator

	output = np.empty(len(initial_state))

	sol = odeint(ODE_function, initial_state, t, args = (output, matrix, k, flux, 
														cross_section))
	if ravel is True:
		sol = np.ravel(sol)

	return sol

def residual_ODE(P, initial, t, y, matrix, numba = True):
	''' Computes residual from data to be fitted and fit '''

	y_fit = ODE_matrix_fit_func(P, initial, t, matrix, numba = numba)
	res = np.sum((y - y_fit)**2)

	return res / 1000.

def ODE_fitting(data, t, reaction_string, sampling, numba = True):

	if numba is True:
		matrix, k_number, reactant_number = reaction_string_to_numba_matrix(reaction_string)
	else:
		matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)

	data_ravel = np.ravel(data)

	bounds_arr = []   # bounds for rate constants in optimization

	for i in range(k_number):
		bounds_arr.append([0.0, np.inf])
	bounds_tuple = tuple(bounds_arr)

	results = []

	for _ in range(sampling):
		guess = np.random.rand(k_number)/10000.
		p = minimize(fun=residual_ODE, x0=guess, args=(data[0], t, data_ravel, matrix, numba), bounds=bounds_tuple, method = 'L-BFGS-B')
		results.append([p.fun, p.x])

	results.sort(key=lambda x: x[0])

	p_ode = results[0][1]

	print("Residual Error:", results[0][0])
	print("Rate Constants:", p_ode)

	return p_ode, matrix

def ODE_conc_profile(p, reaction_string, initial_state, t, plotting = False):

	matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)

	y = ODE_matrix_fit_func(p, initial_state, t, matrix, ravel = False)

	if plotting == True:
		plot_func(y, t, '', show_labels = True, markersize = 10, transpose = True)

	return y

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

	def plot_results(self, labels, ax = plt):

		plot_func(self.data_norm, self.t, '.', show_labels = True, labels = labels, ax = ax, markersize = 10, transpose = True)
		plot_func(self.y_solved, self.t_full, '', ax = ax, markersize = 10, transpose = True)

		if ax != plt:
			ax.set_xlabel('Time / h')
			ax.set_ylabel('Normalized Concentration / a.u.')
			ax.legend()
			ax.set_ylim(-0.05,1)
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	def ODE_plot_generic(self, p, reaction_string, initial_state, t):

		matrix, k_number, reactant_number = reaction_string_to_matrix(reaction_string)

		y = ODE_matrix_fit_func(p, initial_state, t, matrix, ravel = False)

		plot_func(y, t, '', show_labels = True, markersize = 10, transpose = True)

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
	#main()
	plt.show()