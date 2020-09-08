import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.optimize import minimize
from numpy.linalg import inv
import find_nearest as fn
import Liquid_Phase_O2_Analysis as lp
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']


def expo(a, x):

	return a[0] * np.exp(a[1]*x)

def poly_expo(a, x):

	return expo(a[:2], x) + lp.poly(a[2:], x)

def gaussian(a, x):

	return a[0] * np.exp(-((x-a[1])**2)/(a[2]**2))

def voigt(a, x):

	return a[0] * lorentzian(a[1:], x) + (1 - a[0]) * gaussian(a[1:], x)

def voigt_mult(a, x):
	'''Sum of multiple Voigt functions with number of Voigt functions determined by length of a. Length of
	a has to be a multiple of 4.'''

	a_mat = np.reshape(a, (-1, 4))
	y = voigt(a_mat[0], x)

	for i in a_mat[1:,]:
		y = y + voigt(i, x)

	return y

def spline_fitting(x, y, x_full, y_full):
	'''Spline fitting using spline interpolation implemented in SciPy.'''

	tck = interpolate.splrep(x, y, s = 0)
	y_fit = interpolate.splev(x_full, tck, der = 0)
	y_corrected = y_full - y_fit

	return y_fit, y_corrected

def polynomial_regression(x, y, x_full, y_full, order):
	'''Polynomial regression function, where x and y are used for fitting and x_full and y_full represent full
	data range.'''

	X_raw = np.tile(x, (order+1, 1))
	powers = np.arange(0, order+1)

	X = np.power(X_raw.T, powers)

	coef = np.dot(np.dot(inv(np.dot(X.T, X)), X.T), y)

	y_fit = lp.poly(coef, x_full)
	y_corrected = y_full - y_fit
	
	return y_fit, y_corrected

def baseline_fitting_generic(data, feature_start, feature_end, function, p_guess = None, order_poly = 1):
	'''Baseline Correction using either spline_fitting, polynomial_regression (using order_poly) 
	or any other defined function via least_squares (using p_guess).

	se_array provides information on pre-edge peak region beginning, end and signal beginning and end (details see prep_data_fitting function)
	function determines the function used for baseline correction.

	p_guess is used for any functions fitted using least_squares. length of p_guess determines polynomial order or number of voigt peaks, if either of those two functions are used.
	Order_poly is used for polynomial regression.'''

	idx = fn.find_nearest(data, (feature_start, feature_end))

	pre_feature = data[:idx[0]]
	post_feature = data[idx[1]:]
	baseline = np.r_[pre_feature, post_feature]

	if function == spline_fitting:
		y_baseline, y_corrected = spline_fitting(baseline[:,0], baseline[:,1], data[:,0], data[:,1])

	elif function == polynomial_regression:
		y_baseline, y_corrected = polynomial_regression(baseline[:,0], baseline[:,1], data[:,0], data[:,1], order_poly)

	else:
		bounds = np.zeros((2, len(p_guess)))
		bounds[0] = -np.inf
		bounds[1] = np.inf

		if function == voigt:
			bounds[0] = 0.0
			bounds[1][0] = 1.0

		p = least_squares(fun=lp.residual_generic, x0=p_guess, args=(baseline[:,0], baseline[:,1], function), bounds = bounds)
		p_solved = p.x

		y_baseline = function(p_solved, data[:,0])
		y_corrected = data[:,1] - y_baseline

	data_corr = np.c_[data[:,0], y_corrected]

	return data_corr, y_baseline

def integrated_rate_law(k, t):
	'''Integrated rate law for concentration of B in reaction scheme A -> B -> C'''

	return k[2] * (k[0]/(k[1] - k[0])) * (np.exp(-k[0]*t) - np.exp(-k[1]*t))

def residual_rate_law(k, t, y):

	y_fit = integrated_rate_law(k, t)
	res = np.sum((y - y_fit)**2)

	return res

def fit_follow_up(data_raw, start, end, y_multiplier):

	intv = fn.find_nearest(data_raw, (start, end))

	data = data_raw[intv[0]:intv[1]]
	data_x = (data[:,0] - data[0][0])/3600
	data_y = (data[:,1] - data[0][1])*y_multiplier

	#guess = np.ones(3)
	guess = np.array([1., 800., 0.1])

	k_solved = minimize(fun=residual_rate_law, x0 = guess, args=(data_x, data_y), method='Nelder-Mead')
	print(k_solved.success)
	print(k_solved.x)

	k_sol = k_solved.x

	plt.plot(data_x, integrated_rate_law(k_sol, data_x))
	plt.plot(data_x, data_y)

class Experiment:

	def __init__(self, name, feature_start, feature_end, name_b = None, offset = None):

		self.data = lp.import_txt('../Experimental_Data/Gas_Phase_O2_Data/%s' % name, channel = 2, sensors = 1)
		self.feature_start = feature_start
		self.feature_end = feature_end

		if name_b is not None:
			data_b = lp.import_txt('../Experimental_Data/Gas_Phase_O2_Data/%s' % name_b, channel = 2, sensors = 1)
			final_x = self.data[:,0][-1]
			data_b[:,0] = data_b[:,0] + final_x + offset

			self.data = np.r_[self.data, data_b]
			self.feature_start = self.feature_start + final_x + offset

	def fit_baseline(self, function, order_poly = 3, p_guess = None):

		self.data_corr, self.baseline = baseline_fitting_generic(self.data, self.feature_start, self.feature_end, function, order_poly = order_poly, p_guess = None)
	
	def smooth_baseline_corrected_data(self, window_length, poly_order):

		self.y_smooth = savgol_filter(self.data_corr[:,1], window_length, poly_order)

	def fit_follow_up_kinetics(self, y_multiplier = 10000):

		fit_follow_up(self.data_corr, self.feature_start, self.feature_end, y_multiplier)

	def plot_raw_data(self):

		plt.plot(self.data[:,0], self.data[:,1], linewidth = 0.7, color = 'black')

	def plot_baseline_corrected_data(self, only_corrected = True, smoothed = False, offset_correction = False, ax = plt, width = 0.0003):

		if smoothed == True:
			ax.plot(self.data_corr[:,0], self.y_smooth, color = 'green', linewidth = 1.5, label = 'Data')
		else:
			ax.plot(self.data_corr[:,0], self.data_corr[:,1], color = 'black', linewidth = 0.7, label = 'Data')

		if offset_correction == True:
			self.data_corr[:,0] = self.data_corr[:,0] - self.feature_start
			ax.plot((0., 0.), (0. - width, 0. + width), color = 'red', label = 'Irradiation Start', linewidth = 2)
		else:
			ax.plot((self.feature_start, self.feature_start), (0. - width, 0. + width), color = 'red', label = 'Irradiation Start', linewidth = 2)

		if only_corrected is not True:
			ax.plot(self.data[:,0], self.data[:,1])
			ax.plot(self.data[:,0], self.baseline)

		if ax != plt:
			ax.grid(color = 'grey', linestyle = '--', linewidth = 0.2)
			ax.set_xlabel('Time / s')
			ax.set_ylabel(r'$O_{2}$ / Vol%')
			ax.legend()


def main():
	js_552 = Experiment('190619_JS_552_2.txt', 16350., 46300.)
	js_555 = Experiment('190723_JS_555_2.txt', 5000., 47000., name_b = '190723_JS_555_3.txt', offset = 60)

	js_552.fit_baseline(polynomial_regression, order_poly = 6)
	js_555.fit_baseline(polynomial_regression, order_poly = 6)
	js_555.smooth_baseline_corrected_data(101, 3)

	fig, ax = plt.subplots()

	js_552.plot_baseline_corrected_data(only_corrected = False, ax = ax, offset_correction = False)

	#js_555.plot_baseline_corrected_data(ax = ax, smoothed = True, offset_correction = True)

	return js_552


if __name__ == '__main__':
	main()
	plt.show()