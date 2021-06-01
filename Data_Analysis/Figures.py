import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
from utility_functions import insert_image
import Liquid_Phase_O2_Analysis as lp
import Free_Energy_Profile as fe
import Reaction_ODE_Fitting as ode
import Photochemical_ODE as photo
import Gas_Phase_O2_Analysis as gas
import OO_Bond_Scans as oo
import Decay_Associated_Spectra as das
import UV_Vis_Spectra as uvvis
import IR_Analysis as ir
import O2_Data_Plotting as o2
import Solar_to_Hydrogen_Efficiency as sth
import Time_Resolved_UV_Vis as mcr

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

def Dual_Irradiation_Errorbars():

	fig, ax, _, _ = lp.main(errorbars = True)

	insert_image('../Computational_Data/Images/B_mono_image.png', 0.07, 0.45, 0.07, ax[1])
	ax[1].annotate(r'[$\bf{B-Mono-Up}$]D$_{0}$', xy = (0.09, 0.24), xycoords = 'axes fraction')

	return fig

def Dual_Irradiation_f_singlet(spectrum_for_fit = 'f_singlet', legend_loc = 'upper right', split_label_lines = False):

	theoretical_spectra = uvvis.import_theoretical_spectra(return_f_triplet = True)
	spectra_parameters = uvvis.import_plotting_parameters()
	spectrum_for_fit = spectrum_for_fit

	dual_irradiation = lp.main(return_dual_irradiation = True)

	fig, ax = plt.subplots(figsize = (5,4))
	fig.subplots_adjust(left = 0.2, bottom = 0.15, top = 0.95)

	dual_irradiation.fit_theoretical_spectrum_to_synergy(theoretical_spectra[spectrum_for_fit], 'lr', 'lr', 3400., 1000.)
	dual_irradiation.plot_theoretical_fit(ax, spectra_parameters[spectrum_for_fit]['Name'], split_label_lines = split_label_lines)

	ax.legend(loc = legend_loc)

	return fig, ax

def Dual_Irradiation_f_triplet():

	fig, ax = Dual_Irradiation_f_singlet(spectrum_for_fit = 'f_triplet')

	return fig 

def Dual_Irradiation_g_mono():

	fig, ax = Dual_Irradiation_f_singlet(spectrum_for_fit = 'g_mono', legend_loc = 'lower left')

	insert_image('../Computational_Data/Images/B_mono_image.png', 0.04, 0.425, 0.07, ax)
	ax.annotate(r'[$\bf{B-Mono}$]D$_{0}$', xy = (0.09, 0.24), xycoords = 'axes fraction')

	return fig

def Dual_Irradiation_g_trans_mono_three_trans_triplet():

	theoretical_spectra = uvvis.import_theoretical_spectra(return_f_triplet = True)
	spectra_parameters = uvvis.import_plotting_parameters()

	dual_irradiation = lp.main(return_dual_irradiation = True)

	fig, ax = plt.subplots(1,2, figsize = (9,4))
	fig.subplots_adjust(left = 0.1, right = 0.98, bottom = 0.15, top = 0.9, wspace = 0.3)

	spectrum_for_fit = 'g_trans_mono'

	dual_irradiation.fit_theoretical_spectrum_to_synergy(theoretical_spectra[spectrum_for_fit], 'lr', 'lr', 3400., 1000.)
	dual_irradiation.plot_theoretical_fit(ax[0], spectra_parameters[spectrum_for_fit]['Name'], split_label_lines = True)

	spectrum_for_fit = 'three_trans_triplet'

	dual_irradiation.fit_theoretical_spectrum_to_synergy(theoretical_spectra[spectrum_for_fit], 'lr', 'lr', 3400., 1000.)
	dual_irradiation.plot_theoretical_fit(ax[1], spectra_parameters[spectrum_for_fit]['Name'], split_label_lines = True)	

	ax[0].legend(loc = 'upper right')
	ax[1].legend(loc = 'upper right')

	ax[0].text(-0.12, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')

	insert_image('../Computational_Data/Images/B_trans_mono_image.png', 0.03, 0.225, 0.06, ax[0])
	ax[0].annotate(r'[$\bf{B-Mono-Trans}$]D$_{0}$', xy = (0.02, 0.06), xycoords = 'axes fraction')

	insert_image('../Computational_Data/Images/B_trans_mono_image.png', 0.03, 0.225, 0.06, ax[1])
	ax[1].annotate(r'[$\bf{B-Mono-Trans-Up}$]D$_{0}$', xy = (0.02, 0.06), xycoords = 'axes fraction')

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
	ax = profile.final_ax

	return fig

def Free_Energy_Profile_with_ChemDraw():

	profile = fe.Free_Energy_Profile('../Computational_Data/Overview_Energies/DFT_Energies.xlsx')
	fig = profile.final_figure
	ax = profile.final_ax

	ax.text(-0.035, -0.1, 'B', transform=ax.transAxes, size = 20, weight='bold')
	insert_image('../Figures/ChemDraw/Mechanism_no_ox_labels.png', 0.5, -0.27, 0.15, ax)

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

	ode_model, _, _ = photo.main()
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

	fig, ax = das.main()

	insert_image('../Computational_Data/Images/B_trans_image.png', 0.45, 0.19, 0.058, ax)
	ax.annotate(r'[$\bf{B-Trans}$]T$_{0}$', xy = (0.53, 0.025), xycoords = 'axes fraction')

	return fig

def Decay_Associated_Spectra_f_triplet():

	fig, ax = das.secondary()

	return fig

def Decay_Associated_Spectra_g_r2_r2():

	fig, ax = das.tertiary()

	insert_image('../Computational_Data/Images/B_image.png', 0.5, 0.19, 0.052, ax)
	ax.annotate(r'[$\bf{B}$]T$_{0}$', xy = (0.47, 0.025), xycoords = 'axes fraction')

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

	fig, ax = uvvis.secondary()

	insert_image('../Computational_Data/Images/A_image.png', 0.75, 0.5, 0.055, ax)
	ax.annotate(r'[$\bf{A}$]S$_{0}$', xy = (0.72, 0.29), xycoords = 'axes fraction')

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

def STH_Two_Photon():

	fig = sth.secondary()

	return fig

def STH_Three_Photon():

	fig = sth.tertiary()

	return fig

def MCR_Heatmap():

	fig = mcr.MCR_heatmap()

	return fig

def MCR_Results():

	fig = mcr.MCR_results()

	return fig

def UVVis_Data():

	fig = mcr.UVVis_data()

	return fig

def ODE_Fit_to_MCR():

	fig = mcr.ODE_fit_to_MCR()

	return fig

def Intensity_Photo_ODE():

	intensity, actinometry = lp.main(return_intensity = True)
	ode_model, _, _ = photo.main()

	fig, ax = plt.subplots(1, 2, figsize = (9.2,4))
	fig.subplots_adjust(wspace = 0.25, bottom = 0.15, left = 0.09, right = 0.82) 

	intensity.plot_results('fit', False, ax[0], photon_flux = actinometry.p_scaled[0])
	intensity.plot_results('average', False, ax[0], photon_flux = actinometry.p_scaled[0])
	
	ode_model.plot_ODE_fit(ax = ax[1])

	ax[0].text(-0.12, 1.05, 'A', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.05, 'B', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[1].text(1.05, 0.42, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	insert_image('../Figures/ChemDraw/Photochemical_Kinetic_Scheme.png', 1.3, 0.17, 0.22, ax[1])

	return fig

def Bimolecular_Photo_ODE():

	fig, ax = photo.Bimolecular()

	insert_image('../Figures/ChemDraw/Bimolecular_Reaction_Parallel.png', -0.35, 1.26, 0.28, ax[1])

	ax[1].text(-0.87, 1.40, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Three_Photon_Irreversible():

	fig, ax = photo.Three_photon_mechanism_irreversible()

	insert_image('../Figures/ChemDraw/Sequential_Three_Photon.png', -0.35, 1.3, 0.28, ax[1])

	ax[1].text(-0.97, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Three_Photon_Irreversible_Differential_Evolution():

	fig, ax = photo.Three_photon_mechanism_irreversible(differential_evolution_result = True)

	insert_image('../Figures/ChemDraw/Sequential_Three_Photon.png', -0.35, 1.3, 0.28, ax[1])

	ax[1].text(-0.97, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Three_Photon_Reversible():

	fig, ax = photo.Three_photon_mechanism_reversible()

	insert_image('../Figures/ChemDraw/Sequential_Three_Photon_Reversible.png', -0.35, 1.25, 0.28, ax[1])

	ax[1].text(-0.97, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Three_Photon_Reversible_Cubic():

	fig, ax = photo.Three_photon_mechanism_reversible_cubic()

	insert_image('../Figures/ChemDraw/Sequential_Three_Photon_Reversible.png', -0.35, 1.25, 0.28, ax[1])

	ax[1].text(-0.97, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Two_Trans_UV_Vis_Milstein():

	fig = uvvis.two_trans_uv_vis()

	return fig

def Intensity_Linear_Fit():

	fig = lp.flux_rate_relationship_linear_fit()

	return fig

def Simultanous_Two_Photon_Bimolecular():

	fig, ax = photo.simultanous_two_photon_bimolecular()

	insert_image('../Figures/ChemDraw/Bimolecular_Reaction_Two_Photon.png', -0.35, 1.25, 0.28, ax[1])

	ax[1].text(-0.87, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Simultanous_Two_Photon_Bimolecular_Quartic():

	fig, ax = photo.simultanous_two_photon_bimolecular_quartic()

	insert_image('../Figures/ChemDraw/Bimolecular_Reaction_Two_Photon.png', -0.35, 1.25, 0.28, ax[1])

	ax[1].text(-0.87, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

def Simultanous_Two_Photon_Bimolecular_Irreversible():

	fig, ax = photo.simultanous_two_photon_bimolecular_irreversible()

	insert_image('../Figures/ChemDraw/Bimolecular_Reaction_Two_Photon_Irreversible_Thermal.png', -0.35, 1.25, 0.28, ax[1])

	ax[1].text(-0.87, 1.33, 'A', transform=ax[1].transAxes, size = 20, weight='bold')
	ax[0].text(-0.12, 1.02, 'B', transform=ax[0].transAxes, size = 20, weight='bold')
	ax[1].text(-0.12, 1.02, 'C', transform=ax[1].transAxes, size = 20, weight='bold')

	return fig

class Figure:

	def __init__(self, function):
		self.function = function
		self.name = function.__name__
		self.fig = function()

	def show(self):
		plt.show()

	def save(self):
		self.fig.savefig('../Figures/%s.pdf' % self.name, transparent = True, dpi = 700)


figure = Figure(Conical_Intersections)
figure.show()
#figure.save()

