# Overview

This repository contains original data and analysis code for "Two-photon, visible light water splitting at a molecular ruthenium complex".

The publication (open access) is available [here](https://doi.org/10.1039/D1EE01053K) and you can directly access the PDF file [here](https://github.com/jschneidewind/Water-Splitting/blob/master/Documents/Paper_d1ee01053k.pdf).

`Computational_Data` contains computational data for all computed intermediates (typically .log files and visualizations)

`Data_Analysis` contains the different Python programs used for data analysis. Note that these programs use various files from `Computational_Data` as well as `Experimental_Data`, hence these directory must be present to run most of them.

`Documents` contains formatted documents such as the supporting information.

`Experimental_Data` contains various experimental data, especially spectroscopy and kinetic data as well Excel sheets used for small calculations on experimental data.

`Figures` contains figures that appear in either the manuscript or supporting information (typically as .pdf files)

Subdirectories contain seperate README.md files which detail their contents.

We actively encourage others to work on this interesting new avenue in water splitting and are happy to answer questions or comments and are open to cooperations in this field.

# Dependencies

Python code (Python >3.7) requires the following libraries: numba, numpy, scipy, matplotlib, pandas, pathlib, csv and pymcr.

# Installation

To run data analysis code yourself it is easiest to clone the entire repository:

```bash
git clone https://github.com/jschneidewind/Water-Splitting
```


# License

Copyright (c) Jacob Schneidewind

All software is licensed under a MIT license (see `LICENSE` file).

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

All other files and their contents are licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. (see `LICENSE-CC-BY`)

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
