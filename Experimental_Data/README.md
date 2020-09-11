# General

This directory contains experimental data.

`Absorbance_Fluorescence_Data` contains UV/Vis and fluoresence spectra for complex **1**.

`EPR_Data` contains EPR data for complex **1**.

`Gas_Phase_O2_Data` contains gas phase O<sub>2</sub>  measurements for the irradiation of **1**.

`IR_Data` contains IR data for complex **1**.

`Liquid_PHase_O2_Data` contains liquid phase O<sub>2</sub> measurements for the irradiaiton of **1** in kinetic experiments. It contains an Excel sheet with metadata for all experiments (`Liquid_Phase_O2_Experiments_Metadata.xlsx`). Metadata is to be interpreted in the following way:

1. `type` indicates the dataset, which an experiment is part of (`intensity`, `temperature`, `kie` or `dual`)

2. `offset` is the start time for irradiation (in seconds), to which an `offset_correction` might be applied

3. `feature_end` is the end time of irradiation relative to the start time (in seconds) 

4. `wavelength` indicates for irradiation wavelength. For `intensity`, `temperature` and `kie` experiments a wavelength of `500` indicates that the Hg light source (320 - 500 nm) was used. For `dual` experiments, `wavelength` indicates the cut-on wavelength for the used longpass filter

5. `power` indicates the intensity setting of the Hg light source for `intensity` experiments. For `dual` experiments, a power of `0.01` indicates that only the longpass filtered QTH light source was used. A power of `0.02` indicates that only the Hg light source (320 - 400 nm) was used. A power of `0.03` indicates that both light sources were used

6. `active` indicates whether the experiment is used for data analysis

7. `start` indicates an alternate start time for O<sub>2</sub> baseline determination

`NMR_Data` contains NMR data (in `Spectra_Export`, .mnova files and .pdf visualizations) and Excel sheets with numerical data for kinetic and dual irradiation NMR experiments.

`Photographs` contains photographs of experimental set-ups.

`Chemical_Actinometry.xlsx` contains chemical actinometry data.

`Dual_Irradiation_Experimental_Predicted_Rates.xlsx` contains a calculation comparing experimental excess initial rates from dual irradiation experiments with predictions from the photochemical kinetic model.

`ODE_Kinetic_Fit.xlsx` contains optimized parameters for the kinetic model fitted to kinetic NMR data.

`Photochemical_Kinetic_Fit.xlsx` contains optimized parameters for the photochemical kinetic model.

`QTH_Measured_Power.xlsx` contains experimentally measured power of QTH light source with different longpass filters.

`Rate_Comparison_O2_Hydride_Formation.xlsx` contains a calculation comparing the rates of **2-cis** and O<sub>2</sub> formation.