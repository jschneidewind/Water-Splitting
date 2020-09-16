# Overview

This repository contains original data and analysis code for "Two-photon water splitting at a molecular ruthenium complex".

`Computational_Data` contains computational data for all computed intermediates (typically .log files and visualizations)

`Data_Analysis` contains the different Python programs used for data analysis. Note that these programs use various files from `Computational_Data` as well as `Experimental_Data`, hence these directory must be present to run most of them.

`Documents` contains formatted documents such as the supporting information.

`Experimental_Data` contains various experimental data, especially spectroscopy and kinetic data as well Excel sheets used for small calculations on experimental data.

`Figures` contains figures that appear in either the manuscript or supporting information (typically as .pdf files)

Subdirectories contain seperate README.md files which detail their contents.

We actively encourage others to work on this interesting new avenue in water splitting and are happy to answer questions or comments and are open to cooperations in this field.

# Dependencies

Python code (Python >3.7) requires the following libraries: numpy, scipy, matplotlib, pandas, pathlib, math, io, os, re and csv.

# Installation

To run data analysis code yourself it is easiest to clone the entire repository:

```bash
git clone https://github.com/jschneidewind/Water-Splitting
```


