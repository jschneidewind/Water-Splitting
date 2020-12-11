import numpy as np
from scipy import constants as con
from scipy.optimize import minimize
from scipy.optimize import dual_annealing
from scipy.optimize import differential_evolution
from scipy.optimize import shgo
import matplotlib.pyplot as plt
import find_nearest as fn
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

def Boltzmann_distribution(energy_diff, T):
	'''Energy diff in Joule'''

	population_ratio = 1./np.exp(energy_diff/(con.k*T))

	#return population_ratio / (1 + population_ratio)
	return population_ratio

def spectra_integration(data, start, end):
	'''Integration of provided 2D data array from start to end values'''

	idx = fn.find_nearest(data, (start, end))
	data_intv = data[idx[0]:idx[1]]
	integral = np.trapz(data_intv[:,1], data_intv[:,0])

	return integral

def nm(value):
	'''Converts nm to J'''
	return ((con.h*con.c)/(value/1e9))

def eV(value):
	'''Converts eV to J'''
	return value*1.602176e-19

def kcalmol(value):
	'''Converts kcalmol to J'''
	return (value * 4186.798188)/con.Avogadro

def Jmol(value):
	'''Convert Jmol to J'''
	return value / con.Avogadro

class Energy:

	def __init__(self, value, unit):
		'''Input value in either nm, eV or kcalmol
		Available units: J, eV, nm, kcal/mol and J/mol'''
		self.unit = unit.__name__
		self.value = value

		self.J = unit(self.value)
		self.eV = self.convert_J_to_eV()
		self.nm = self.convert_J_to_nm()
		self.kcalmol = self.convert_J_to_kcalmol()
		self.Jmol = self.convert_J_to_Jmol()

	def convert_J_to_eV(self):
		return self.J/1.602176e-19

	def convert_J_to_nm(self):
		return ((con.h*con.c)/self.J)*1e9

	def convert_J_to_kcalmol(self):
		return (self.J * con.Avogadro)/4186.798188

	def convert_J_to_Jmol(self):
		return self.J * con.Avogadro

class Efficiency:

	def __init__(self, gibbs_energy, excess_energy):
		'''Initiating AM 1.5 G spectra (irradiance, flux and total irradiace) as well as Gibbs energy for considered
		reaction and excess energy (energy that is lost due to free energy and kinetic losses)
		gibbs_energy and excess_energy are energy class instances'''

		self.am15g_full = self.import_am15g_data()
		self.am15g_irradiance = np.c_[self.am15g_full[:,0], self.am15g_full[:,2]]
		self.am15g_flux = self.convert_spectrum_to_photon_count(self.am15g_irradiance)	
		self.total_irradiance = spectra_integration(self.am15g_irradiance, 0., 5000.)
		self.gibbs_energy = gibbs_energy
		self.excess_energy = excess_energy

	def import_am15g_data(self):
		'''AM1.5G ASTM G-173 reference data from https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html'''
		data = np.loadtxt('../Experimental_Data/ASTMG173.csv', skiprows = 2, delimiter = ',')
		return data

	def convert_spectrum_to_photon_count(self, spectrum):
		'''Input spectrum in nm, W m^-2 nm^-1
		returns photon count in photons m^-2 s^-1 nm^-1'''

		photon_energies = Energy(spectrum[:,0], nm)
		photon_count = spectrum[:,1]/photon_energies.J

		return np.c_[spectrum[:,0], photon_count]

	def STH_Efficiency(self,reactions_per_second, delta_G, irradiance, area):
		'''Solar-to-hydrogen efficiency, reactions_per_second in mol/s, delta_G in J/mol,
		irradiance in W/m^2, area in m^2

		If delta_G is 237 kJ/mol, reactions_per_second refers to H2 mol/s
		If delta_G is 474 kJ/mol,reactions_per_second refers to 2 H2 mol/s (input value is halve)
		Arbitrary delta_G value can be sued to account for partial water splitting reaction'''

		return (reactions_per_second * delta_G)/(irradiance * area)

	def calc_efficiency(self, wavelengths, number_photons, negative = False, return_reactions = False):
		'''Consecutive absorption of photons by a number of intermediates. Number of intermediates determined by
		dimension of wavelengths array.
		Intermediate that absorbs the lowest number of photons is the bottleneck and determines the number of 
		water splitting reactions per s.
		number_photons is an array of equal length as wavelengths, specifying the number of photons absorbed at corresponding
		wavelength
		Based on excess_energy and provided wavelengths, the amount of available energy (variable "energy") is calculated.
		if the available energy is lower than gibbs_energy, the population of sufficiently energtic particles is calculated
		using a two state Boltzmann model and this population is multiplied with min(photons_per_second_list) to obtain reactions_per_second'''
	
		energy = Energy(Energy(wavelengths, nm).eV - self.excess_energy.eV, eV)
		energy_difference = Energy(self.gibbs_energy.eV - np.sum(number_photons * energy.eV), eV)

		if energy_difference.eV > 0:
			pop = Boltzmann_distribution(energy_difference.J, 298.15)
		else:
			pop = 1.

		photons_per_second_list = []

		for counter, (wavelength, photons) in enumerate(zip(wavelengths, number_photons)):
			if counter == 0:
				photons_per_second = spectra_integration(self.am15g_flux, 0., wavelength) / con.Avogadro
			else:
				photons_per_second = spectra_integration(self.am15g_flux, wavelengths[counter-1], wavelength) / con.Avogadro
			
			photons_per_second_list.append(photons_per_second/photons)

		reactions_per_second = pop * min(photons_per_second_list)

		sth = self.STH_Efficiency(reactions_per_second, self.gibbs_energy.Jmol, self.total_irradiance, 1.)

		if return_reactions is True:
			return reactions_per_second
		elif negative is True:
			return -sth
		else:
			return sth

	def plot_1d(self, wavelengths, number_photons):

		sth_values = []

		for i in wavelengths:
			sth = self.calc_efficiency(np.array([i]), number_photons)
			sth_values.append(sth)

		sth_values = np.asarray(sth_values)

		fig, ax = plt.subplots()
		ax.plot(wavelengths, sth_values)

		return fig, ax

	def plot_2d_colormesh(self, wavelengths_a, wavelengths_b, number_photons, print_maximum = False):

		x, y = np.meshgrid(wavelengths_a, wavelengths_b)

		sth_values = []

		for i in zip(np.ravel(x), np.ravel(y)):
			sth = self.calc_efficiency(np.array([i[0], i[1]]), number_photons)
			sth_values.append(sth)

		sth_values = np.asarray(sth_values)

		if print_maximum is True:
			idx = np.argmax(sth_values)
			print('Maximum STH:', sth_values[idx])
			print('Optimal wavelengths:', np.ravel(x)[idx], np.ravel(y)[idx])

		sth_values = np.reshape(sth_values, x.shape)

		levels = MaxNLocator(nbins=200).tick_values(np.amin(sth_values), np.amax(sth_values))
		cmap = plt.get_cmap('inferno')
		norm = BoundaryNorm(levels, ncolors=cmap.N, clip = True)

		fig, ax = plt.subplots()

		im_a = ax.pcolormesh(x, y, sth_values, cmap=cmap, norm=norm)
		im_a.set_edgecolor('face')
		colorbar = fig.colorbar(im_a, ax=ax)

		colorbar.set_label('STH / %')
		ax.set_xlabel(r'Longest absorption wavelength [$\bf{A}$] / nm')
		ax.set_ylabel(r'Longest absorption wavelength [$\bf{B}$] / nm')

		return fig, ax

	def plot_solar_irradiance(self, start, end, ax = plt):

		idx_plot = fn.find_nearest(self.am15g_flux, (start, end))
		data = self.am15g_flux[idx_plot[0]:idx_plot[1]]

		ax.plot(data[:,0], data[:,1], color = 'black', linewidth = 1., label = 'AM 1.5 G')

		idx_segments = np.asarray(fn.find_nearest(self.am15g_flux, self.optimal_wavelengths))
		full_idx_segments = np.insert(idx_segments, 0, 200)
		midpoints = (full_idx_segments[1:] + full_idx_segments[:-1]) / 2
		midpoints = midpoints.astype(int)

		offsets = np.array([-1, 0, 1]) * 80.

		segments = np.split(self.am15g_flux, idx_segments)

		cmap = plt.get_cmap('inferno')
		color_list = cmap(np.linspace(0, 1, 5))

		for counter, (i, idx, idx_mid, wavelength, offset) in enumerate(zip(segments[:-1], idx_segments, midpoints, self.optimal_wavelengths, offsets)):
			
			photon_flux = spectra_integration(i, 0., 5000.)
			i = np.vstack((i, self.am15g_flux[idx]))

			ax.fill_between(i[:,0], i[:,1], color = color_list[counter+1])
			ax.annotate(r'$\lambda_{%s}$ =' % (counter + 1) + '\n' + ' %.0f nm' % wavelength + '\nFlux $(m^{-2} s^{-1})$ =\n' + r'${0:s}$'.format(as_si(photon_flux, 2)), 
					(self.am15g_flux[:,0][idx_mid] + offset, 5e18), ha = 'center', color = color_list[counter])

		ax.annotate('Maximum STH:\n' + r'$\bf{%.2f}$' % (self.optimal_efficiency * 100) + r'$\bf{\%}$',
				    (0.81, 0.73), color = 'black', ha = 'center', xycoords = 'figure fraction')

		if ax != plt:
			ax.set_xlabel('Wavelength / nm')
			ax.set_ylabel(r'Spectral photon flux / $m^{-2} nm^{-1} s^{-1}$')
			ax.legend()
			ax.set_ylim(0, 6.2e18)

	def locate_optimum(self, p_guess, number_photons, method = 'Nelder-Mead', print_output = False):

		p = minimize(fun=self.calc_efficiency, x0=p_guess, args=(number_photons, True), method = method)

		self.optimal_wavelengths = p.x
		self.optimal_efficiency = -p.fun
		self.number_photons = number_photons

		if print_output is True:
			print('Photon sequence:', number_photons, 'Wavelengths (nm):', self.optimal_wavelengths, 'STH (%):', '%.2f' % (100.*self.optimal_efficiency))

	def global_optimization(self, number_photons, method = 'differential evolution', print_output = False):

		bounds = np.array([300., 3000.])
		bounds = np.tile(bounds, (len(number_photons),1))

		if method == 'dual annealing':
			p = dual_annealing(func=self.calc_efficiency, bounds = bounds, args = (number_photons, True), maxiter = 1000, local_search_options = {'method': 'Nelder-Mead'})
		
		elif method == 'shgo':
			p = shgo(func=self.calc_efficiency, bounds = bounds, args =(number_photons, True))

		else:
			p = differential_evolution(func=self.calc_efficiency, bounds = bounds, args = (number_photons, True))

		self.optimal_wavelengths = p.x
		self.optimal_efficiency = -p.fun
		self.number_photons = number_photons

		if print_output is True:
			print('Photon sequence:', number_photons, 'Wavelengths (nm):', np.around(self.optimal_wavelengths, 0) , 'STH (%):', '%.2f' % (100.*self.optimal_efficiency))

	def monte_carlo_optimization(self, number_photons, samples, print_results = False):

		results = []

		for _ in range(0, samples):

			p_guess = np.random.uniform(300., 3000., number_photons.shape)
			p_guess = p_guess[np.argsort(p_guess)]
			print(p_guess)
			self.locate_optimum(p_guess, number_photons)
			sth = self.optimal_efficiency

			results.append(sth)

		results = np.asarray(results)
		results = results[np.argsort(results)]

		if print_results is True:
			print(results)

	def thermal_contribution(self, complete_gibbs_energy, print_output = False):
		
		reactions_per_second = self.calc_efficiency(self.optimal_wavelengths, self.number_photons, return_reactions = True)
		thermal_energy = complete_gibbs_energy.Jmol - self.gibbs_energy.Jmol
		
		sth = self.STH_Efficiency(reactions_per_second, thermal_energy, self.total_irradiance, 1.)

		if print_output is True:
			print('Photon sequence:', self.number_photons, 'Wavelengths (nm):', self.optimal_wavelengths, 'STH thermal contribution (%):', '%.2f' % (100.*sth))
			print('Total STH (%):', '%.2f' % (100. * (sth + self.optimal_efficiency)))

def main():
	#gibbs_energy = Energy(82., kcalmol)
	gibbs_energy = Energy(4*1.229, eV)
	excess_energy = Energy(17.5, kcalmol)
	#excess_energy = Energy(1., eV)

	efficiency = Efficiency(gibbs_energy, excess_energy)

	efficiency.locate_optimum(p_guess = np.array([450.]), number_photons = np.array([2.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450.]), number_photons = np.array([4.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550.]), number_photons = np.array([1., 1.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550.]), number_photons = np.array([2., 2.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550.]), number_photons = np.array([4., 4.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550., 650.]), number_photons = np.array([1., 1., 1.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550., 650.]), number_photons = np.array([2., 2., 2.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550., 650.]), number_photons = np.array([1., 1., 2.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550., 650., 750.]), number_photons = np.array([1., 1., 1., 1.]), print_output = True)
	efficiency.locate_optimum(p_guess = np.array([450., 550., 650., 750.])*2., number_photons = np.array([2., 2., 2., 2.]), print_output = True)
	#efficiency.locate_optimum(p_guess = np.array([450., 550., 650., 1750.]), number_photons = np.array([2., 2., 2., 2.]), print_output = True)

def secondary():
	gibbs_energy = Energy(82.7, kcalmol)
	excess_energy = Energy(17.5, kcalmol)
	#excess_energy = Energy(1., eV)

	efficiency = Efficiency(gibbs_energy, excess_energy)

	wavelengths_a = np.linspace(300., 500., 100)
	wavelengths_b = np.linspace(300., 1000., 100)

	#print(spectra_integration(efficiency.am15g_flux, 0., 455.)/10000.)

	fig, ax = efficiency.plot_2d_colormesh(wavelengths_a, wavelengths_b, number_photons = np.array([1., 1.]))

	efficiency.locate_optimum(p_guess = np.array([450., 550.]), number_photons = np.array([1., 1.]))


	efficiency.thermal_contribution(Energy(117.11, kcalmol), print_output = True)

	ax.annotate('Maximum STH: %.2f' % (efficiency.optimal_efficiency * 100) + '%' + '\n%.0f nm' % efficiency.optimal_wavelengths[0] +
		', %.0f nm' % efficiency.optimal_wavelengths[1], (0.6, 0.8), color = 'white', ha = 'center', xycoords = 'figure fraction')

	return fig

def tertiary():
	gibbs_energy = Energy(4*1.229, eV)
	excess_energy = Energy(17.5, kcalmol)

	efficiency = Efficiency(gibbs_energy, excess_energy)

	efficiency.locate_optimum(p_guess = np.array([450., 550., 650.]), number_photons = np.array([1., 1., 2.]), print_output = True)

	fig, ax = plt.subplots()

	efficiency.plot_solar_irradiance(200., 1000., ax = ax)

	return fig

def quaternary():

	gibbs_energy = Energy(4*1.229, eV)
	#excess_energy = Energy(17.5, kcalmol)
	excess_energy = Energy(1., eV)

	efficiency = Efficiency(gibbs_energy, excess_energy)
	efficiency.monte_carlo_optimization(np.array([1., 1., 1., 1., 1.])*1, 100, print_results = True)

def quinary():

	method = 'differential evolution'

	gibbs_energy = Energy(4*1.229, eV)
	excess_energy = Energy(1., eV)
	#excess_energy = Energy(17.5, kcalmol)

	efficiency = Efficiency(gibbs_energy, excess_energy)

	efficiency.global_optimization(np.array([2.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([4.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([1., 1.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([2., 2.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([4., 4.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([1., 1., 1.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([2., 2., 2.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([1., 1., 2.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([1., 1., 1., 1.]), method = method, print_output = True)
	efficiency.global_optimization(np.array([2., 2., 2., 2.]), method = method, print_output = True)


if __name__ == '__main__':
	#main()
	secondary()
	#tertiary()
	#quaternary()
	#quinary()
	plt.show()
