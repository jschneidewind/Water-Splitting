import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import Liquid_Phase_O2_Analysis as lp
import Free_Energy_Profile as fe
import ODE_Fit_to_Data as ode
import Photochemical_Reaction_ODE_Sim as photo
import Gas_Phase_O2_Analysis as gas
import OO_Bond_Scans as oo
import Decay_Associated_Spectra as das
import UV_Vis_Spectra as uvvis
import IR_Analysis as ir
import O2_Data_Plotting as o2

def Intensity_Dual_Irradiation():
	fig, ax, _, _ = lp.main(dual_and_intensity = True)

	fig.subplots_adjust(wspace = 0.3, bottom = 0.15)

	ax[0].text(-0.07, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.07, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Intensity():

	_, _, fig, ax = lp.main()

	return fig

def Dual_Irradiation():

	fig, ax, _, _ = lp.main()

	return fig

def Dual_Irradiation_f_singlet(spectrum_for_fit = 'f_singlet', legend_loc = 'upper right'):

	theoretical_spectra = uvvis.import_theoretical_spectra(return_f_triplet = True)
	spectra_parameters = uvvis.import_plotting_parameters()
	spectrum_for_fit = spectrum_for_fit

	dual_irradiation = lp.main(return_dual_irradiation = True)

	fig, ax = plt.subplots(figsize = (5,4))
	fig.subplots_adjust(left = 0.2, bottom = 0.15, top = 0.95)

	dual_irradiation.fit_theoretical_spectrum_to_synergy(theoretical_spectra[spectrum_for_fit], 'lr', 'lr', 3400., 1000.)
	dual_irradiation.plot_theoretical_fit(ax, spectra_parameters[spectrum_for_fit]['Name'])

	ax.legend(loc = legend_loc)

	return fig

def Dual_Irradiation_f_triplet():

	fig = Dual_Irradiation_f_singlet(spectrum_for_fit = 'f_triplet')

	return fig 

def Dual_Irradiation_g_mono():

	fig = Dual_Irradiation_f_singlet(spectrum_for_fit = 'g_mono', legend_loc = 'lower left')

	return fig

def H2O2_Disproportionation():
	fig, ax = lp.secondary()

	fig.subplots_adjust(left = 0.08, right = 0.92)

	ax[0].text(-0.07, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.07, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Chemical_Actinometry():

	actinometry = lp.Chemical_Actinometry('../Experimental_Data/Chemical_Actinometry.xlsx', '../Experimental_Data/20120613_Lumatec2_Spektren.txt')
	fig, ax = plt.subplots()
	actinometry.plot_linear_regression(ax = ax)

	return fig

def Blackbody_Fit():

	dual_irradiation = lp.main(return_dual_irradiation = True)
	fig, ax = plt.subplots()
	dual_irradiation.abs_spectrum_dual.check_blackbody_fit(ax = ax)

	return fig

def Free_Energy_Profile():

	profile = fe.Free_Energy_Profile('../Computational_Data/Overview_Energies/DFT_Energies.xlsx')
	fig = profile.final_figure

	return fig

def Conical_Intersections():

	conical = fe.Conical_Intersections('../Computational_Data/Overview_Energies/CASSCF_Energies.xlsx')
	fig = conical.final_figure

	return fig

def Conical_Intersections_PCM():

	conical = fe.Conical_Intersections('../Computational_Data/Overview_Energies/CASSCF_Energies.xlsx', pcm = True)
	fig = conical.final_figure

	return fig

def ODE_fit():

	kinetic_model, labels = ode.main()
	fig, ax = plt.subplots()
	kinetic_model.plot_results(labels, ax = ax)

	return fig

def Photo_ODE_Model():

	ode_model = photo.main()
	fig, ax = plt.subplots(1,2, figsize = (11.5, 4.5))
	fig.subplots_adjust(wspace = 0.6)

	ode_model.plot_ODE_fit(ax = ax[0])
	ode_model.plot_tau_initial_rate(ax = ax[1])

	ax[0].text(-0.07, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.07, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Gas_Phase_O2():

	gas_phase = gas.main()
	fig, ax = plt.subplots()
	fig.subplots_adjust(left = 0.15)
	gas_phase.plot_baseline_corrected_data(only_corrected = True, ax = ax)

	return fig

def OO_Scans():

	fig = oo.main()

	return fig

def Decay_Associated_Spectra():

	fig = das.main()

	return fig

def Decay_Associated_Spectra_f_triplet():

	fig = das.secondary()

	return fig

def Theoretical_UV_Vis_Spectra():

	fig = uvvis.main()

	return fig

def DE_Scan():

	fig = oo.secondary()

	return fig

def IR_Spectrum():

	fig = ir.main()

	return fig

def EPR_Spectrum():

	fig = ir.secondary()

	return fig

def UV_Vis_Fluorescence():

	fig = uvvis.secondary()

	return fig

def Lumatec_Spectra():

	fig = uvvis.tertiary()

	return fig

def Liquid_Phase_Raw():

	fig = o2.main()

	return fig

def Liquid_Phase_Analysis():

	fig = lp.tertiary()

	return fig


class Figure:

	def __init__(self, function):
		self.function = function
		self.name = function.__name__
		self.fig = function()

	def show(self):
		plt.show()

	def save(self):
		self.fig.savefig('../Figures/%s.pdf' % self.name, transparent = True, dpi = 500)


figure = Figure(Free_Energy_Profile)
figure.show()
#figure.save()

