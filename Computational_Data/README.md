# General

Directories `A - F` contain corresponding data for each intermediate/transition state/scan. Within each directory there might be multiple subdirectories for each model chemistry:

`Full_Model` for dimeric model chemistry, `Mono` for monomeric model chemistry and `Me_Mono` for monomeric model chemistry with simplified ligand model. If only one model chemistry was used, the directory will simply contain the corresponding files.

For each structure there a typically files corresponding to two or three different calculations (DFT calculations):

1. Geometry optimization and frequency calculation, designated by `-GPOF` file name (or `-GPF` if file contains only the frequency calculation)

2. Single point energy calculation including SMD solvation, designated by `-SMD` file name

3. TD-DFT calculation including SMD solvation, designated by `-TD-SMD` file name. When a TD-DFT calculation was performed, there is typically also a `_uvvis.txt` file in the same directory containing data for the computed UV/Vis spectrum.

For CASSCF calculations, optimizations are designated by `-Opt` file name, while calculations including solvation are designated by `-PCM` file name.

Furthermore, for each intermediate there is an `.xyz` file containing coordinates, a `.png` file containing a visualization and a `.pov` file which can be used to render the visualization.

`Images` contains simplified images used for visualization.

`Misc` contains various other calculations, such as O-O bond scans, results for the **Oxo Dimer** and NMR reference data.

`Overview_Energies` contains Excel sheets in which important energy values are collected and compared. Importantly, `DFT_Energies.xlsx` is used to generate the free energy profile and `CASSCF_Energies.xlsx` is used to generate the reaction path from [**B**]D<sub>0</sub> to [**C**]D<sub>0</sub>.

