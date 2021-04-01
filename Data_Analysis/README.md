# General

This directory contains all Python programs used for data analysis. Most programs depend on each other and import code from other modules, so that for easy use, all programs should be present to run data analysis.

Overview of contents:

Name | Contents
--- | ---
`Absorption_Spectrum_to_Dual_Irradiation.py` | Code to convert theoretical UV/Vis Spectra to predicted dual irradiation behaviour
`Controlled_Average.py` | Module to average sub-sets of data so that intervals between to data sets are identical
`Decay_Associated_Spectra.py` | Code to visualize experimental decay associated spectra (DAS) and to compute theoretical DAS from UV/Vis spectra
`Figures.py` | Code to generate final figures
`find_nearest.py` | Module to find index of value closest to provided one in an array
`Free_Energy_Profile.py` | Code for visualization of free energy profile and conical intersections plot
`Gas_Phase_O2_Analysis.py` | Code for gas phase O<sub>2</sub> analysis
`IR_Analysis.py` | Visualization of experimental and theoretical IR spectra as well as EPR data
`Liquid_Phase_O2_Analysis.py` | Main code for kinetic analysis. Contains code for analysis of intensity and dual irradiation data sets, KIE and temperature dependence, chemical actinometry and H<sub>2</sub>O<sub>2</sub> disproportionation.
`O2_Data_Plotting.py` | Short script to visualize any experimental O<sub>2</sub> measurement file
`OO_Bond_Scans.py` | Visualization of O-O bond scan results as well as Scan-[**DE**]T<sub>0</sub> results
`Photochemical_ODE.py` | Photochemical kinetic model and fitting of model to experimental data (uses ODE code from `Reaction_ODE_Fitting.py`)
`Reaction_ODE_Fitting.py` | Code for constructing and solving ODE system of chemical reaction network,  fitting of kinetic model to NMR data
`Solar_to_Hydrogen_Efficiency.py` | Model for calculation of maximum theoretical solar-to-hydrogen efficiencies
`Time_Resolved_UV_Vis.py` | Analysis of time resolved UV-Vis data using mutlivariate curve resolution (MCR)
`utility_functions.py` | General purpose functions (e.g. plotting, image insertion etc.) used in multiple modules
`UV_Vis_Spectra.py` | Import, processing and visualization of theoretical and experimental UV/Vis as well as light source spectra