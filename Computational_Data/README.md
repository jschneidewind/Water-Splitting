# General

Directories `A - F` contain corresponding data for each intermediate/transition state/scan. Within each directory there might be multiple subdirectories for each model chemistry:

`Full_Model` for dimeric model chemistry, `Mono` for monomeric model chemistry and `Me_Mono` for monomeric model chemistry with simplified ligand model. If only one model chemistry was used, the directory will simply contain the corresponding files.

For each structure there a typically files corresponding to two or three different calculations (DFT calculations):

1. Geometry optimization and frequency calculation, designated by `-GPOF` file name (or `-GPF` if file contains only the frequency calculation)

2. Single point energy calculation including SMD solvation, designated by `-SMD` file name

3. TD-DFT calculation including SMD solvation, designated by `-TD-SMD` file name. When a TD-DFT calculation was performed, there is typically also a `_uvvis.txt` file in the same directory containing data for the computed UV/Vis spectrum.

For DFT calculations, `-GP` indicates the use of the standard functional/basis set combinations (PBE0 + cc-pVDZ or cc-pVTZ). For CASSCF calculations, optimizations are designated by `-Opt` file name, while calculations including solvation are designated by `-PCM` file name.

Furthermore, for each intermediate there is an `.xyz` file containing coordinates, a `.png` file containing a visualization and a `.pov` file which can be used to render the visualization.

`Images` contains simplified images used for visualization.

`Misc` contains various other calculations, such as O-O bond scans, results for the **Oxo Dimer** and NMR reference data.

`Overview_Energies` contains Excel sheets in which important energy values are collected and compared. Importantly, `DFT_Energies.xlsx` is used to generate the free energy profile and `CASSCF_Energies.xlsx` is used to generate the reaction path from [**B**]D<sub>0</sub> to [**C**]D<sub>0</sub>.

# Naming

The names used in manuscript and supporting information differ from those of the original computational files provided herein. The correspondence is as follows:

Original Name | Final Name
--- | ----
F-R2-R2-Down-Disp | [**A**]S<sub>0</sub>
3-Full-Disp | [**A-Mono**]S<sub>0</sub>
G-Cis-R2-R2-Down-T- | [**B**]T<sub>0</sub>
GMe-Cis-Mono-D-0 | [**B-Mono**]D<sub>0</sub> (Me model)
G-Cis-Mono-D-0 | [**B-Mono**]D<sub>0</sub>
3-HCis-D-0 | [**B-Mono-Up**]D<sub>0</sub>
G-R2-R2-Down | [**B-Trans**]T<sub>0</sub>
G-Mono-D-0 | [**B-Mono-Trans**]D<sub>0</sub>
3-HTrans-D-0 | [**B-Mono-Trans-Up**]D<sub>0</sub>
GMe-Cis-CI12 | [**BC**] D<sub>1</sub>/D<sub>0</sub> MECI
GMe-Cis-CI23 | [**BC**] D<sub>2</sub>/D<sub>1</sub> MECI
H-Trans-R2-R2-Down-RuMove-T | [**C**]T<sub>0</sub>
GMe-Cis-OOScan6-D0 *and* HMe-Trans-Mono-D-0 | [**C-Mono**]D<sub>0</sub> (Me model)
H-Trans-Mono-D-0-180 | [**C-Mono**]D<sub>0</sub>
TSHJ-Trans-R2-R2-Down-T | TS-[**CD**]T<sub>0</sub>
J-Trans-R2-R2-Down-TSGuess-T | [**D**]T<sub>0</sub> [**F**]D<sub>0</sub> (Dimer model)
J-Trans-Mono-H2O-T0 | [**D**]T<sub>0</sub>
J-Trans-Mono-H2O-T0-RuO-Scan | Scan-[**DE**]T<sub>0</sub>
0-Full-Planar-H2O | [**E**]S<sub>0</sub>
2-Full-Cis-A | [**F-Cis**]S<sub>0</sub>
2-Full-Cis-B | [**F-Cis-Up**]S<sub>0</sub>
2-Full-TPy-A | [**F**]S<sub>0</sub>
2-Full-Tpy-B | [**F-Up**]S<sub>0</sub>
2-Full-Trans-A | [**F-Trans**]S<sub>0</sub>
2-Full-Trans-B | [**F-Trans-Up**]S<sub>0</sub>
GAlt-Mono | [**B-Alt**]
3-Dimer | [**Oxo Dimer**]S<sub>0</sub>






