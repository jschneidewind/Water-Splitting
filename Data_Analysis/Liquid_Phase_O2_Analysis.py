import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import io as io
import find_nearest as fn
from UV_Vis_Spectra import Constructed_UV_Vis_Spectrum, import_theoretical_spectra, import_plotting_parameters
from Absorption_Spectrum_to_Dual_Irradiation import Absorption_Spectrum_Dual_Irradiation
from import_lamp_spectra import import_lamp_spectra
from scipy.optimize import least_squares
from scipy.stats import linregress
from scipy import constants as con
from pathlib import Path
import pandas as pd
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def autolabel(rects, ax = plt):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}%'.format(int(height)),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize = 12)

def import_txt(name, channel, sensors):

	raw_data = np.genfromtxt(name, skip_header=14+(2*sensors-2), encoding='ISO8859', usecols = (0,1,2,3,4,5,6,7,8,9,10,11))
	data = raw_data[:,[2, 2+channel, 6+channel]]
	
	return data

def poly(a, x):

    y = a[0] * x**0

    for i in range(1, len(a)):
        y += a[i] * x**i

    return y

def integrated_rate_law_expanded(k, initial, t):
	'''Integrated form of d[O2]/dt = k0 - k1 * [O2] for [O2] (t = 0) = initial'''

	return initial * np.exp(-k[1] * t) - k[0] * np.exp(-k[1] * t) * (1/k[1]) + k[0]/k[1]

def integrated_rate_law_expanded_single(k0, initial, t, k1):
	'''Integrated form of d[O2]/dt = k0 - k1 * [O2] for [O2] (t = 0) = initial'''

	return initial * np.exp(-k1 * t) - k0 * np.exp(-k1 * t) * (1/k1) + k0/k1

def integrated_rate_law_shift(p, initial, t, k1 = None):

	idx = fn.find_nearest(t, p[-1])[0]

	t_base = t[:idx+1]
	t_feature = t[idx+1:]

	y_base = t_base * 0

	if k1 is None:
		y_feature = integrated_rate_law_expanded(p[:-1], initial, (t_feature - p[-1]))

	else:
		y_feature = integrated_rate_law_expanded_single(p[:-1], initial, (t_feature - p[-1]), k1)

	return np.r_[y_base, y_feature]

def baseline_feature_function(p, t, order_poly, k1 = None):

	#initial = poly(p[:order_poly+1], p[-1])

	return poly(p[:order_poly+1], t) + integrated_rate_law_shift(p[order_poly+1:], 0, t, k1)

def feature_function(p, t, initial, k1 = None):

	if k1 is None:
		return integrated_rate_law_expanded(p, initial, t)
	else:
		return integrated_rate_law_expanded_single(p, initial, t, k1)

def square_law(p, x):

	return p * x**2

def linear_model(p, x):

	return p * x

def linear_regression(p, x):

	return p[0] * x + p[1]

def power_law(p, x, power):

	return p * x**power

def power_model(a, x):

	return a[0] * x**a[1]

def residual_baseline_feature(p, t, order_poly, y, k1 = None):

	y_fit = baseline_feature_function(p, t, order_poly, k1)
	res = y - y_fit

	return res

def residual_feature(p, t, initial, y, k1 = None):

	y_fit = feature_function(p, t, initial, k1)
	res = y - y_fit

	return res

def residual_generic(p, x, y, function):

	y_fit = function(p, x)
	res = y - y_fit

	return res

def residual_generic_power(p, x, y, power):

	y_fit = power_law(p, x, power)
	res = y - y_fit

	return res

def fit_generic(x, y, function, p_guess = None):

	if p_guess is None:
		if function == linear_regression or function == power_model:
			p_guess = np.ones(2)

		else:
			p_guess = np.ones(1)

	p = least_squares(fun=residual_generic, x0=p_guess, args=(x, y, function))

	x_full = np.linspace(np.amin(x), np.amax(x), 1000.)

	y_solved_full = function(p.x, x_full)
	y_solved = function(p.x, x)

	return x_full, y_solved_full, p.x

def linear_regression_correction(data, column):

	slope, intercept, r_value, p_value, std_err = linregress(data[:,0], data[:,column])
	data_x = np.unique(data[:,0])
	data_lin_reg_values = slope * data_x + intercept

	return data_x, data_lin_reg_values, np.mean(data[:,column])

def consecutive(data, stepsize=1):

	return np.split(data, np.where(np.diff(data[:,0]) != 0)[0]+1)	
				
def sort_average_by_parameter(parameter, rates, are_lists = True, add_zeros = True):

	if are_lists is True:
		parameter = np.asarray(parameter)
		rates = np.asarray(rates)

	data = np.c_[parameter, rates]
	idx_sort = np.argsort(data[:,0], kind = 'mergesort')
	data = data[idx_sort]

	if add_zeros is True:
		data = np.concatenate((np.zeros(3)[None,:], data), axis = 0)

	consec = consecutive(data)

	average_rates = []

	for i in consec:
		average_rates.append(np.average(i, axis = 0))

	return np.asarray(average_rates), data

def average_data(exps):

	power_levels = []
	data = []
	time = []

	for i in exps:
		if exps[i].active == True:
			exps[i].liquid_phase_O2_analysis_feature(only_data_processing = True)
			power_levels.append(exps[i].power)
			data.append(exps[i].feature_O2)
			time.append(exps[i].feature_t)

	power_levels = np.asarray(power_levels)
	data = np.asarray(data)
	time = np.asarray(time)

	idx_sort = np.argsort(power_levels, kind = 'mergesort')

	power_levels = power_levels[idx_sort]	

	data = data[idx_sort]
	time = time[idx_sort]

	idx_split = np.where(np.diff(power_levels) != 0)[0]+1

	data_split = np.split(data, idx_split)
	time_split = np.split(time, idx_split)

	result = []

	for i in data_split:

		data_int = np.zeros(399)

		for j in i:
			if len(j) < 399:
				add = np.random.normal(np.average(j[-20:]), np.std(j[-20:]), 399-len(j))
				j = np.append(j, add)
			data_int += j

		data_avg = data_int / len(i)
		result.append(data_avg)

	return np.unique(power_levels), np.asarray(result), time[-1]

def subtraction_analysis(exps):

	a = 'js_607_7'
	b = 'js_607_8'
	c = 'js_607_9'

	print(exps[a].rate)
	print(exps[b].rate)
	print(exps[c].rate)

	sub = exps[c].data - exps[a].data - exps[c].data

	plt.plot(exps[a].time, sub)

def convert_xlsx_to_experiments(xlsx):

	df = pd.read_excel(xlsx, index_col = 0, dtype = object)

	fill = {'active': True, 'analysis': False, 'color': 'black', 'plotting': False, 'start': 0., 'offset_correction': 0., 'k1': 'None'}

	df = df.fillna(fill)
	dic = df.to_dict('index')

	experiments = {'dual': {}, 'intensity': {}, 'kie': {}, 'temperature': {}}

	for i in dic:

		if dic[i]['k1'] == 'None':
			dic[i]['k1'] = None

		exp = Experiment(name = dic[i]['name'], offset = dic[i]['offset'], feature_end = dic[i]['feature_end'], wavelength = dic[i]['wavelength'],
						 power = dic[i]['power'], color = dic[i]['color'], start = dic[i]['start'], analysis = dic[i]['analysis'], offset_correction = dic[i]['offset_correction'],
						 plotting = dic[i]['plotting'], k1 = dic[i]['k1'], active = dic[i]['active'], exp_type = dic[i]['type'])

		experiments[dic[i]['type']][i] = exp

	return experiments

class Experiment:

	def __init__(self, name, offset, feature_end, wavelength, power, color, start = 0., analysis = True, offset_correction = 0., plotting = False, k1 = None, active = True, exp_type = None):
		self.name = name
		self.offset = offset
		self.start = start
		self.feature_end = feature_end
		self.wavelength = wavelength
		self.power = power
		self.color = color
		self.analysis = analysis
		self.active = active
		self.offset_correction = offset_correction
		self.correct_offset = self.offset + self.offset_correction
		self.k1 = k1
		self.plotting = plotting
		self.exp_type = exp_type

		p = Path('../Experimental_Data/Liquid_Phase_O2_Data/')
		file = p / self.name
		self.raw_data = import_txt(file, channel = 1, sensors = 1)

	def liquid_phase_O2_analysis_combined(self, order_poly, start = None, plotting = None, k1 = None):

		if k1 is None:
			k1 = self.k1
		if start is None:
			start = self.start
		if plotting is None:
			plotting = self.plotting

		data = self.raw_data
		idx = fn.find_nearest(data, (start, self.offset + self.feature_end))
		data = data[idx[0]:idx[1]]

		if k1 is None:
			p_guess = np.ones(order_poly + 4)
		else:
			p_guess = np.ones(order_poly + 3)

		p_guess[-1] = self.offset

		p = least_squares(fun=residual_baseline_feature, x0=p_guess, args=(data[:,0], order_poly, data[:,1], k1), method='lm')

		if k1 is None:
			k_solved = p.x[-3:-1]
		else:
			k_solved = np.array([p.x[-2], k1])

		if plotting == True:
			y_fit = baseline_feature_function(p.x, data[:,0], order_poly, k1)	
			plt.plot(data[:,0], data[:,1])
			plt.plot(data[:,0], y_fit)
			plt.plot(data[:,0], y_fit - data[:,1])

		self.rate_combined = k_solved[0]
		self.rate_combined_full = k_solved

	def liquid_phase_O2_analysis_feature(self, offset = None, start = None, plotting = None, k1 = None, only_data_processing = False, ax = plt):

		if offset is None:
			offset = self.correct_offset
		if start is None:
			start = self.start
		if k1 is None:
			k1 = self.k1
		if plotting is None:
			plotting = self.plotting

		data = self.raw_data
		idx = fn.find_nearest(data, (start, offset, offset + self.feature_end))

		pre_feature_average = np.average(data[:,1][idx[0]:idx[1]])

		feature = data[idx[1]+1:idx[2]]

		feature_t = feature[:,0] - offset
		feature_O2 = feature[:,1] - pre_feature_average

		if only_data_processing == True:
			self.feature_t = feature_t
			self.feature_O2 = feature_O2

		else:
			pre_feature_average = 0

			if k1 is None:
				p_guess = np.ones(2)/1000.
			else:
				p_guess = np.ones(1)/1000.

			p = least_squares(fun=residual_feature, x0=p_guess, args=(feature_t, pre_feature_average, feature_O2, k1))

			if k1 is None:
				k_solved = p.x
			else:
				k_solved = np.array([p.x[0], k1])

			residual = residual_feature(p.x, feature_t, pre_feature_average, feature_O2, k1)

			if plotting == True:
				y_solved = feature_function(p.x, feature_t, pre_feature_average, k1)
				ax.plot(feature_t, feature_O2, 'o-', markersize = 2., linewidth = 0.2, color = self.color, label = 'Data')
				ax.plot(feature_t, y_solved, color = self.color, label = 'Fit')
				ax.plot(feature_t, (feature_O2 - y_solved), color = self.color, linewidth = 0.3, label = 'Residual')

				if ax is not plt:
					ax.set_xlabel('Time / s')
					ax.set_ylabel(r'$O_2$ / $\mu$mol $l^{-1}$')
					ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

				ax.legend()

			self.residual = np.sum(residual[:50]**2)
			self.rate_feature = k_solved[0]
			self.rate_feature_full = k_solved

	def determine_offset_correction(self):

		offsets = np.arange(self.offset, self.offset + 40., 1.)
		residuals = []

		for i in offsets:
			self.liquid_phase_O2_analysis_feature(offset = i, k1 = self.k1)
			residuals.append(self.residual)

		residuals = np.asarray(residuals)
		sort_idx = np.argsort(residuals)

		plt.plot(offsets-self.offset, residuals, 'o-')
		plt.plot(offsets[sort_idx[:5]]-self.offset, residuals[sort_idx[:5]], '.', color = 'red', markersize = 15)

	def plot_raw_data(self):

		plt.plot(self.raw_data[:,0], self.raw_data[:,1], linewidth = 0.7, color = self.color)

class Intensity_Analysis:

	def __init__(self, analysis_function, experiments, order_poly = 0):
		self.analysis_function = analysis_function
		self.exps = experiments
		self.order_poly = order_poly
		self.intensity_analysis()

	def intensity_analysis(self):
		'''Intensity analysis of liquid phase O2 data. Analysis is done using either combined (analysis_function = combined) or feature (analysis_function =
		feature) liquid_phase_O2_analysis function. For combined analysis function, order_poly has to be specified.
		Experiments are provided as a dictionary'''

		power_levels = []
		rates_full = []

		for i in self.exps:
			if self.exps[i].active == True:
				power_levels.append(self.exps[i].power)

				if self.analysis_function == 'combined':
					self.exps[i].liquid_phase_O2_analysis_combined(self.order_poly)
					rates_full.append(self.exps[i].rate_combined_full)

				if self.analysis_function == 'feature':
					self.exps[i].liquid_phase_O2_analysis_feature()
					rates_full.append(self.exps[i].rate_feature_full)

		average_rates, full_rates = sort_average_by_parameter(power_levels, rates_full)

		k1_dict = dict(zip(average_rates[:,0], average_rates[:,2]))

		rates_full_fixed_k1 = [] 

		for i in self.exps:
			if self.exps[i].active == True:

				if self.analysis_function == 'combined':
					self.exps[i].liquid_phase_O2_analysis_combined(self.order_poly, k1 = k1_dict[self.exps[i].power])
					rates_full_fixed_k1.append(self.exps[i].rate_combined_full)

				if self.analysis_function == 'feature':
					self.exps[i].liquid_phase_O2_analysis_feature(k1 = k1_dict[self.exps[i].power])
					rates_full_fixed_k1.append(self.exps[i].rate_feature_full)

		average_rates_fixed_k1, full_rates_fixed_k1 = sort_average_by_parameter(power_levels, rates_full_fixed_k1)

		self.intensity_data_avg = average_rates
		self.intensity_data_full = full_rates

		self.intensity_data_avg_fixed_k1 = average_rates_fixed_k1
		self.intensity_data_full_fixed_k1 = full_rates_fixed_k1

	def return_data(self, data_type, fixed_k1):

		if data_type == 'average':
			if fixed_k1 == False:
				return self.intensity_data_avg
			else:
				return self.intensity_data_avg_fixed_k1

		elif data_type == 'full':
			if fixed_k1 == False:
				return self.intensity_data_full
			else:
				return self.intensity_data_full_fixed_k1

		else:
			pass

	def fit_intensity_data(self, function, data_type, fixed_k1):

		data = self.return_data(data_type, fixed_k1)

		x = data[:,0]
		y = data[:,1]

		x_full, y_solved, p = fit_generic(x, y, function)

		print(p)

		self.fit = np.c_[x_full, y_solved]
		self.fit_function = function

	def log_linear_analysis(self, data_type, fixed_k1):

		data = self.return_data(data_type, fixed_k1)
		data = data[1:,]

		x = data[:,0]
		y = np.log10(data[:,1])

		x_full, y_solved, p = fit_generic(x, y, linear_regression)

		self.log_linear_fit = np.c_[x_full, y_solved]

	def plot_results(self, data_type, fixed_k1, ax = plt, photon_flux = None):
		'''Plotting function for intensity analysis results. data_type specifies if full or averaged data is plotted. fixed_k1
		(True or False) specifies if data is plotted, which was obtained by fixing k1 or not'''

		if data_type == 'average':
			if fixed_k1 == False:
				ax.plot(self.intensity_data_avg[:,0], self.intensity_data_avg[:,1], '.', markersize = 25, color = 'black', label = 'Datapoints')
			else:
				ax.plot(self.intensity_data_avg_fixed_k1[:,0], self.intensity_data_avg[:,1], '.', markersize = 10, color = 'darkorange')

		if data_type == 'full':
			if fixed_k1 == False:
				ax.plot(self.intensity_data_full[:,0], self.intensity_data_full[:,1], '.', markersize = 5, color = 'lightgreen')
			else:
				ax.plot(self.intensity_data_full_fixed_k1[:,0], self.intensity_data_full_fixed_k1[:,1], '.', markersize = 5, color = 'orange')

		if data_type == 'fit':
			if self.fit_function == square_law:
				ax.plot(self.fit[:,0], self.fit[:,1], color = 'green', label = r'Fit using $f(x) = ax^{2}$', linewidth = 2)
			else:
				ax.plot(self.fit[:,0], self.fit[:,1], color = 'black')

		if data_type == 'log_linear':
			ax.plot(self.log_linear_fit[::,0], self.log_linear_fit[:,1], color = 'black')

		if ax != plt:

			ax.set_xlabel(r'Photon Flux (320 - 500 nm) / $ s^{-1}$')
			ax.set_ylabel(r'Initial Rate of $O_2$ Formation / $\mu mol s^{-1} l^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()

			f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
			g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
			fmt = mticker.FuncFormatter(g)

			multiplier = '{:.2f}'.format(photon_flux/1e18)
			#ax.set_xlabel(r'Photon Flux (320 - 500 nm) / $ s^{-1} \cdot$' + '{:.2e}'.format(photon_flux))
	
			ax.text(0.85, -0.12, r'$\times$' + multiplier + r'$\cdot$' + '{}'.format(fmt(1e18)), transform=ax.transAxes)

class Dual_Irradiation_Analysis:

	def __init__(self, analysis_function, experiments, order_poly = 0):
		self.analysis_function = analysis_function
		self.exps = experiments
		self.order_poly = order_poly
		self.dual_irradiation_analysis()

	def import_QTH_data(self, name):

		df = pd.read_excel(name)
		data = df.to_numpy()

		self.filter_wavelengths = data[:,0]
		self.measured_power = data[:,1]

	def dual_irradiation_analysis(self):
		'''Analysis of dual irradiation data. Three types of experiments are used: QTH irradiation with different long pass
		filters (indicated by single_var and 0.01 power level), Hg irradiation 320 - 400 nm (single_320 and 0.02 power level),
		and dual irradiation (combination of Hg irradiation 320 - 400 nm + variable QTH irradiation, dual, 0.03 power level).

		Analysis is performed using either combined or feature liquid_phase_O2_analysis (for combined analysis, order of poly
		has to be specified).

		First, k0 and k1 values are determined for all experiments. Then, the mean value of all k1 values for each type of experiment
		is determined. Also, a linear regression model is used for single_var and dual experiments to estimate k1 values for each
		wavelength.

		liquid_phase_O2_analysis is again performed, this time fixing k1 to be either the mean value for each type of experiment or fixing
		k1 to be the corresponding lin_reg value at corresponding wavelength, obtaining k0_mean and k0_lr values for each type of experiment.

		The obtained k0 values are then analyzed to determine the synergy present in dual irradiation experiments.
		For this, the k0 values for single_var experiments are estimated at each wavelength using a linear regression model
		or the mean value across all wavelenghts is calculated.

		The corresponding linear regression value at each wavelength, or the mean value, of single_var is subtracted from
		the averaged dual k0 values. Furthermore, the mean k0 value for single_320 experiments is also subtracted.

		Ultimately, the synergy value at each wavelength is obtained in four ways:

		1. k1 linear regression estimation, single_var k0 linear regression estimation
		2. k1 linear regression estimation, sinlge_var k0 mean estimation
		3. k1 mean estimation, single_var k0 linear regression estimation
		4. k1 mean estimation, single_var k0 mean estimation.

		For completeness, synergy is also calculated by subtracting the averaged single_var k0 values from dual k0 (designated as indi)

		'''

		single_var = []
		single_320 = []
		dual = []

		for i in self.exps:
			if self.exps[i].active == True:

				if self.analysis_function == 'combined':
					self.exps[i].liquid_phase_O2_analysis_combined(self.order_poly)
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_combined_full]

				if self.analysis_function == 'feature':
					self.exps[i].liquid_phase_O2_analysis_feature()
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_feature_full]

				if self.exps[i].power == 0.01:
					single_var.append(result)

				if self.exps[i].power == 0.02:
					single_320.append(result)

				if self.exps[i].power == 0.03:
					dual.append(result)

		single_var = np.asarray(single_var)
		single_320 = np.asarray(single_320)
		dual = np.asarray(dual)

		single_var_wavelengths, single_var_lin_reg_values, single_var_mean = linear_regression_correction(single_var, 2)
		single_320_mean = np.mean(single_320[:,2])
		dual_wavelengths, dual_lin_reg_values, dual_mean = linear_regression_correction(dual, 2)

		var_dict = dict(zip(single_var_wavelengths, single_var_lin_reg_values))
		dual_dict = dict(zip(dual_wavelengths, dual_lin_reg_values))
		single_320_dict = {320: single_320_mean}

		k1_values_mean = {0.01: single_var_mean, 0.02: single_320_mean, 0.03: dual_mean}
		k1_values_lin_reg = {0.01: var_dict, 0.02: single_320_dict, 0.03: dual_dict}

		k0_single_var_mean = []  
		k0_single_320_mean = []  
		k0_dual_mean = []

		k0_single_var_lr = []
		k0_dual_lr = []

		mean_arrays = {0.01: k0_single_var_mean, 0.02: k0_single_320_mean, 0.03: k0_dual_mean}
		lr_arrays = {0.01: k0_single_var_lr, 0.02: k0_single_320_mean, 0.03: k0_dual_lr}

		for i in self.exps:
			if self.exps[i].active == True:

				if self.analysis_function == 'combined':
					self.exps[i].liquid_phase_O2_analysis_combined(self.order_poly, k1 = k1_values_mean[self.exps[i].power])
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_combined_full]
					mean_arrays[self.exps[i].power].append(result)

					self.exps[i].liquid_phase_O2_analysis_combined(self.order_poly, k1 = k1_values_lin_reg[self.exps[i].power][self.exps[i].wavelength])
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_combined_full]
					lr_arrays[self.exps[i].power].append(result)

				if self.analysis_function == 'feature':
					self.exps[i].liquid_phase_O2_analysis_feature(k1 = k1_values_mean[self.exps[i].power])
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_feature_full]
					mean_arrays[self.exps[i].power].append(result)

					self.exps[i].liquid_phase_O2_analysis_feature(k1 = k1_values_lin_reg[self.exps[i].power][self.exps[i].wavelength])
					result = np.r_[self.exps[i].wavelength, self.exps[i].rate_feature_full]
					lr_arrays[self.exps[i].power].append(result)

		syn_lr_mean, syn_mean_mean, syn_indi_mean = self.process_k0_values(mean_arrays)
		syn_lr_lr, syn_mean_lr, syn_indi_lr = self.process_k0_values(lr_arrays)

		self.synergy = {'mean': {'lr': syn_lr_mean, 'mean': syn_mean_mean, 'indi': syn_indi_mean},
						'lr': {'lr': syn_lr_lr, 'mean': syn_mean_lr, 'indi': syn_indi_lr}}

	def process_k0_values(self, k0_dict):
		'''Processing dictionary of k0 values and calculating synergy using either linear regression, mean or individual estimation
		for single_var k0 values
		'''

		k0_single_var = np.asarray(k0_dict[0.01])
		k0_single_320 = np.asarray(k0_dict[0.02])
		k0_dual = np.asarray(k0_dict[0.03])

		single_var_wavelengths, k0_single_var_lin_reg_values, k0_single_var_mean = linear_regression_correction(k0_single_var, 1)
		k0_single_320_mean = np.mean(k0_single_320[:,1])

		k0_single_var_avg, k0_single_var_sort = sort_average_by_parameter(k0_single_var[:,0], k0_single_var[:,1], are_lists = False, add_zeros = False)
		k0_dual_avg, k0_dual_sort = sort_average_by_parameter(k0_dual[:,0], k0_dual[:,1], are_lists = False, add_zeros = False)

		synergy_lr = k0_dual_avg[:,1] - k0_single_var_lin_reg_values - k0_single_320_mean
		synergy_mean = k0_dual_avg[:,1] - k0_single_var_mean - k0_single_320_mean
		synergy_individual = k0_dual_avg[:,1] - k0_single_var_avg[:,1] - k0_single_320_mean

		wavelength = k0_dual_avg[:,0]

		return np.c_[wavelength, synergy_lr], np.c_[wavelength, synergy_mean], np.c_[wavelength, synergy_individual]

	def fit_theoretical_spectrum_to_synergy(self, abs_spectrum, k1_correction, k0_correction, lamp_temperature, ir_cutoff):

		self.import_QTH_data('../Experimental_Data/QTH_Measured_Power.xlsx')

		data = self.synergy[k1_correction][k0_correction]

		abs_spectrum_dual = Absorption_Spectrum_Dual_Irradiation(self.filter_wavelengths, self.measured_power, lamp_temperature, ir_cutoff, data[:,0])
		abs_spectrum_dual.convert_absorption_spectrum_to_rate(abs_spectrum)

		self.abs_spectrum_dual = abs_spectrum_dual
		self.exp_rates_fitted = data
		self.theoretical_rates_scaled = abs_spectrum_dual.fit_theoretical_rates_to_experimental(data)

	def plot_theoretical_fit(self, ax = plt, name_fit = 'None'):

		ax.plot(self.theoretical_rates_scaled[:,0], self.theoretical_rates_scaled[:,1], color = 'green', label = r'Predicted for %s (scaled)' % name_fit, linewidth = 2)
		ax.plot(self.exp_rates_fitted[:,0], self.exp_rates_fitted[:,1], '.', color = 'black', markersize = 25, label = 'Datapoints')
		
		#print(self.exp_rates_fitted)

		if ax != plt:
			ax.set_xlabel('Longpass Cut-On Wavelength Second Light Source / nm')
			ax.set_ylabel(r'Excess Initial Rate of $O_2$ Formation / $\mu mol$ $s^{-1}$ $l^{-1}$')
			ax.legend(loc = 'lower left')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

	def plot_synergy_results(self, k1_correction, k0_correction, ax = plt):

		data = self.synergy[k1_correction][k0_correction]

		#print(data)

		ax.plot(data[:,0], data[:,1], '.')

	def plot_hydride_yields(self, ax = plt):

		yields = np.array([[0., 6.], [1., 0.], [2., 10.1]])
		colors = ['darkblue', 'darkred', 'darkgreen']
		labels = ['320 - 400 nm', '> 495 nm', '320 - 400 nm +\n> 495 nm']

		for counter, i in enumerate(yields):
			rect = ax.bar(i[0], i[1], color = colors[counter])
			autolabel(rect, ax = ax)

		if ax != plt:
			ax.set_ylim(0.0, 11.)
			ax.set_ylabel(r'Yield of $\bf{2-cis}$ / %')
			ax.set_xticks(yields[:,0])
			ax.set_xticklabels(labels)

class KIE_Temperature_Analysis:

	def __init__(self, analysis_function, experiments_temperature, experiments_kie, reference_experiments, full_reference, order_poly = 0):
		self.analysis_function = analysis_function
		self.exps_t = experiments_temperature
		self.exps_kie = experiments_kie
		self.exps_ref = reference_experiments
		self.order_poly = order_poly
		self.full_reference = full_reference
		self.kie_temperature_analysis()

	def kie_temperature_analysis(self):
		'''Analysis of H/D KIE and effect of temperature on the reaction (by comparing k_30C/k_20C). Reference are standard reactions
		performed during the same measurement run (Intensity = 1W).

		k1 mean values are determined by optimizing both k0 and k1. Subsequently, a new fit is performed fixing k1 to be
		the mean value for each group of experiments.'''

		rates = {'intensity': [], 'kie': [], 'temperature': []}

		self.determine_rates(self.exps_t, rates)
		self.determine_rates(self.exps_kie, rates)
		self.determine_rates(self.exps_ref, rates)

		k1_mean_values = {}

		for i in rates:
			k1_mean = np.mean(np.asarray(rates[i]), axis = 0)[1]
			k1_mean_values[i] = k1_mean

		rates_k1_fixed = {'intensity': [], 'kie': [], 'temperature': []}

		self.determine_rates(self.exps_t, rates_k1_fixed, k1 = k1_mean_values['temperature'])
		self.determine_rates(self.exps_kie, rates_k1_fixed, k1 = k1_mean_values['kie'])
		self.determine_rates(self.exps_ref, rates_k1_fixed, k1 = k1_mean_values['intensity'])

		k0_mean_values = {}

		for i in rates_k1_fixed:
			k0_mean = np.mean(np.asarray(rates_k1_fixed[i]), axis = 0)[0]
			k0_mean_values[i] = k0_mean

		self.hd_kie = k0_mean_values['intensity'] / k0_mean_values['kie']
		self.k30_k20 = k0_mean_values['temperature'] / k0_mean_values['intensity'] 

	def determine_rates(self, exps, dic, k1 = None):

		for i in exps:

			if exps[i].exp_type == 'temperature' or exps[i].exp_type == 'kie':
				self.call_liquid_phase_O2_analysis(exps[i], dic, k1)

			if exps[i].exp_type == 'intensity':

				if self.full_reference == True:
					if exps[i].power == 1 and exps[i].active == True:
						self.call_liquid_phase_O2_analysis(exps[i], dic, k1)

				else:
					if str.split(exps[i].name, '_')[0] == '200328':
						self.call_liquid_phase_O2_analysis(exps[i], dic, k1)

	def call_liquid_phase_O2_analysis(self, exp, dic, k1):

		if self.analysis_function == 'combined':
			exp.liquid_phase_O2_analysis_combined(self.order_poly, k1 = k1)
			dic[exp.exp_type].append(exp.rate_combined_full)

		if self.analysis_function == 'feature':
			exp.liquid_phase_O2_analysis_feature(k1 = k1)
			dic[exp.exp_type].append(exp.rate_feature_full)

class Chemical_Actinometry:

	def __init__(self, file_name, spectra_name):
		self.file_name = file_name
		self.spectra_name = spectra_name
		self.import_actinometry_data()
		self.spectra_320_400, self.spectra_320_500 = import_lamp_spectra(self.spectra_name)
		self.determine_photon_flux()
		self.linear_regression_flux()

	def import_actinometry_data(self):

		df = pd.read_excel(self.file_name, dtype = object)
		dic = df.to_dict()

		self.power_320_400 = dic['Power (W, LP Standard Head)'][0]
		self.power_320_500 = dic['320 - 500 Power with 320 - 400 Intensity Setting (W, LP Standard Head)'][0]
		self.photon_data = df[['Power (W, LP Standard Head)', 'Photons/s']].to_numpy()

	def determine_photon_flux(self):

		photon_spectra_320_500 = self.convert_spectrum_to_photon_count(self.spectra_320_500)

		int_to_390 = self.spectra_integration(photon_spectra_320_500, 320, 390)
		int_to_500 = self.spectra_integration(photon_spectra_320_500, 320, 500)

		experimental_ratio = self.power_320_400/self.power_320_500
		power_ratio = self.spectra_integration(self.spectra_320_500, 320, 390)/self.spectra_integration(self.spectra_320_500, 320, 500)
		spectral_ratio = int_to_390/int_to_500

		difference = experimental_ratio/power_ratio

		scaling_factor = self.photon_data[:,1][0]/int_to_390
		photon_spectra_320_500_scaled = np.c_[photon_spectra_320_500[:,0], scaling_factor * photon_spectra_320_500[:,1]]

		photon_flux_320_500 = self.spectra_integration(photon_spectra_320_500_scaled, 320, 500)

		self.power_photon_flux_320_500 = np.array([self.power_320_500, photon_flux_320_500])

	def linear_regression_flux(self):

		data = self.photon_data[1:]
		data = np.concatenate((data, np.zeros(2)[None,:]), axis = 0).astype('float')

		x_full, y_fit, p = fit_generic(data[:,0], data[:,1], linear_model, p_guess = 1e17)
		scaling_factor = self.power_photon_flux_320_500[1]/(p*self.power_photon_flux_320_500[0])
		p_scaled = p * scaling_factor

		print(p_scaled)

		self.photon_data_320_500 = data
		self.fit = np.c_[x_full, y_fit]
		self.p = p
		self.p_scaled = p_scaled

	def spectra_integration(self, data, start, end):

		idx = fn.find_nearest(data, (start, end))
		data_intv = data[idx[0]:idx[1]]
		integral = np.trapz(data_intv[:,1], data_intv[:,0])

		return integral

	def convert_spectrum_to_photon_count(self, spectrum):

		photon_energy = self.convert_nm_to_eV(spectrum[:,0])
		photon_count = spectrum[:,1]/photon_energy

		return np.c_[spectrum[:,0], photon_count]

	def convert_nm_to_eV(self, wavelength):
		'''Energy in eV, wavelength in nm'''

		return ((con.h*con.c)/(wavelength/1e9))/(1.602*10**(-19))

	def plot_linear_regression(self, ax = plt):

		ax.plot(self.fit[:,0], self.fit[:,1], linewidth = 3, color = 'green', label = 'Linear Fit')
		ax.plot(self.photon_data_320_500[:,0], self.photon_data_320_500[:,1], '.', markersize = 30, color = 'black', label = 'Datapoints (320 - 500 nm)')
		
		if ax != plt:
			ax.set_xlabel('Measured Power / W')
			ax.set_ylabel(r'Photon Flux / s$^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()

class H2O2_Disproportionation:

	def __init__(self, name, start, feature_start, feature_end):
		self.name = name
		self.start = start
		self.feature_start = feature_start
		self.feature_end = feature_end

		p = Path('../Experimental_Data/Liquid_Phase_O2_Data/')
		file = p / self.name
		self.data = import_txt(file, channel = 1, sensors = 1)

	def disproportionation_analysis(self):

		idx = fn.find_nearest(self.data, (self.start, self.feature_start, self.feature_end))

		self.equilibration = self.data[:idx[0]]
		self.disp_phase = self.data[idx[0]:idx[1]]
		self.feature = self.data[idx[1]:idx[2]]
		self.post_feature = self.data[idx[2]:]

		self.disp_full_t, self.disp_fit, self.disp_p = fit_generic(self.disp_phase[:,0], self.disp_phase[:,1], linear_regression)
		self.post_full_t, self.post_fit, self.post_p = fit_generic(self.post_feature[:,0], self.post_feature[:,1], linear_regression)

	def plot_results(self, ax = plt):

		ax.plot(self.equilibration[:,0], self.equilibration[:,1], label = 'Equilibration', color = 'grey', linewidth = 0.7, alpha = 0.3)
		ax.plot(self.disp_phase[:,0], self.disp_phase[:,1], '.', label = 'Before KI Addition', color = 'orange', markersize = 4)
		ax.plot(self.disp_full_t, self.disp_fit, label = 'Before KI Addition Linear Fit', color = 'orange', linewidth = 2)
		ax.plot(self.feature[:,0], self.feature[:,1], label = 'KI Addition', color = 'red', linewidth = 0.7, alpha = 0.3)
		ax.plot(self.post_feature[:,0], self.post_feature[:,1], '.', label = 'After KI Addition', color = 'green', markersize = 4)
		ax.plot(self.post_full_t, self.post_fit, label = 'After KI Addition Linear Fit', color = 'green', linewidth = 2)

		print(self.disp_p)
		print(self.post_p)

		if ax != plt:
			ax.set_ylim(0.00, 0.08)
			ax.set_xlim(400, 5750)
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_{2}$ / $\mu$mol $l^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.legend()

	def plot_raw_data(self, ax = plt):

		ax.plot(self.data[:,0], self.data[:,1], linewidth = 0.7, color = 'black')

		if ax != plt:
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_{2}$ / $\mu$mol $l^{-1}$')
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)

def main(return_dual_irradiation = False, dual_and_intensity = False):

	excel_exps = convert_xlsx_to_experiments('../Experimental_Data/Liquid_Phase_O2_Data/Liquid_Phase_O2_Experiments_Metadata.xlsx')

	fig, ax = plt.subplots(1, 2, figsize = (9.2,4))
	fig.subplots_adjust(wspace = 0.3, bottom = 0.15, left = 0.07, right = 0.97) 

	actinometry = Chemical_Actinometry('../Experimental_Data/Chemical_Actinometry.xlsx', '../Experimental_Data/20120613_Lumatec2_Spektren.txt')

	# gme_casscf = Constructed_UV_Vis_Spectrum('../Computational_Data/B_Intermediate/Me_Mono/CASSCF/GMe-Cis-D0-VDZ-SA-CAS1311-Pi-RASSI.xlsx', 0.33)
	theoretical_spectra = import_theoretical_spectra()
	spectra_parameters = import_plotting_parameters()
	spectrum_for_fit = 'three_hcis'

	dual_irradiation = Dual_Irradiation_Analysis('feature', excel_exps['dual'], order_poly = 1)
	dual_irradiation.fit_theoretical_spectrum_to_synergy(theoretical_spectra[spectrum_for_fit], 'lr', 'lr', 3400., 1000.)
	dual_irradiation.plot_theoretical_fit(ax[1], spectra_parameters[spectrum_for_fit]['Name'])

	intensity = Intensity_Analysis('feature', excel_exps['intensity'], order_poly = 1)
	intensity.fit_intensity_data(square_law, 'full', False)

	if dual_and_intensity is False:
		dual_irradiation.plot_hydride_yields(ax = ax[0])

	else:
		intensity.plot_results('fit', False, ax[0], photon_flux = actinometry.p_scaled[0])
		intensity.plot_results('average', False, ax[0], photon_flux = actinometry.p_scaled[0])

	ax[0].text(-0.07, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.07, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	fig_b, ax_b = plt.subplots(figsize = (5,4))
	fig_b.subplots_adjust(left = 0.2, bottom = 0.15, top = 0.95)

	intensity.plot_results('fit', False, ax_b, photon_flux = actinometry.p_scaled[0])
	intensity.plot_results('average', False, ax_b, photon_flux = actinometry.p_scaled[0])

	kie_temp = KIE_Temperature_Analysis('feature', excel_exps['temperature'], excel_exps['kie'], excel_exps['intensity'], True, order_poly = 1)
	print(kie_temp.hd_kie)
	print(kie_temp.k30_k20)

	if return_dual_irradiation is True:
		return dual_irradiation

	else:
		return fig, ax, fig_b, ax_b

def secondary():

	fig, ax = plt.subplots(1, 2, figsize = (11,4))

	h2o2_disp = H2O2_Disproportionation('200808_JS_614_1.txt', 700., 3000., 3450.)
	h2o2_disp.disproportionation_analysis()
	h2o2_disp.plot_raw_data(ax[0])
	h2o2_disp.plot_results(ax[1])

	return fig, ax

def tertiary():

	fig, ax = plt.subplots()

	excel_exps = convert_xlsx_to_experiments('../Experimental_Data/Liquid_Phase_O2_Data/Liquid_Phase_O2_Experiments_Metadata.xlsx')
	excel_exps['intensity']['js_600_4'].liquid_phase_O2_analysis_feature(plotting = True, ax = ax)

	return fig


if __name__ == '__main__':
	main()
	#secondary()
	#tertiary()
	plt.show()