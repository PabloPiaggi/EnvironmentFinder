# <img src="https://github.com/PabloPiaggi/EnvironmentFinder/raw/master/logo-notext.png" width="64" valign="middle" alt="EnviornmentFinder"/> Environment finder

A tool for finding and analyzing atomic environments in crystal structures.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4746325.svg)](https://doi.org/10.5281/zenodo.4746325)

## Use the tool online with Binder

<!-- Launch the tool! [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PabloPiaggi/EnvironmentFinder/master?urlpath=apps%2FApp.ipynb) -->
Launch the tool! [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PabloPiaggi/EnvironmentFinder/voila-test?urlpath=voila%2Frender%2FApp.ipynb)

## Run locally

You can run the app locally with conda and the following commands:
```
git clone https://github.com/PabloPiaggi/EnvironmentFinder EnvironmentFinder
cd EnvironmentFinder
conda create --name env_finder python=3.8 anaconda
conda activate env_finder
pip install -r requirements.txt
bash postBuild
jupyter-notebook App.ipynb
```
Once the jupyter notebook opens in your browser, click on "Appmode", and start using the tool.

## Purpose

The purpose of this tool is to find atomic environments in crystal structures.
These environments can then be used to define measures of similarity and collective variables for enhanced sampling simulations.
The output from this tool can be used directly to create reference environments for the [EnvironmentSimilarity](https://www.plumed.org/doc-master/user-doc/html/_e_n_v_i_r_o_n_m_e_n_t_s_i_m_i_l_a_r_i_t_y.html) collective variable in [PLUMED](https://www.plumed.org/doc-master/user-doc/html/index.html).
The tool can also be used to visualize chemical environments around an atom for general purposes.

## Features

* Choose from the examples or upload your own configuration.
* Determine unique environments in crystal structures. Generally, these environments correspond to the *basis* of the crystal structure.
* An algorithm to find unique environments is used.
* Visualize the environments.
* Output the environments in Protein Data Bank (pdb) format.

## Instructions

The app has four tabs:
* Choose configuration
* Define environments
* Analyze environments
* Output environments

that should be used sequentially.

### Choose configuration
The *Choose configuration* tab has a dropdown menu and an *Upload* button.
By default the user can select a configuration from some examples that are useful to play with and get a grasp of the tool.
The *Upload* button should be used to analyze a configuration stored in the user's hard drive.
Once a configuration file has been uploaded the file should appear in the dropdown choices.
Any file format supported by the ASE library can be used.
A complete list can be found [here](https://wiki.fysik.dtu.dk/ase/ase/io/io.html).
ASE determines the file format by its extension so be sure to name the file properly!

### Define environments
The *Define environments* tab is where the user defines the environments.
There are three different ways to define the environments, *Type*, *String*, and *Step* that can be selected with the toggle buttons.
Typically the most useful one is the *Type* option.
In this case one chooses the atom type of the central atom, the atom type of the neighboring atoms, and a cutoff.
These parameters are set using the *Central atoms type:*, *Neighbor atoms type:*, and *Cutoff (Å):* boxes.
The app will find all environments with atoms of the second type around atoms of the first type inside the specified cutoff.

The tool then calculates the unique environments.
For this purpose it compares all environments and keeps only the ones that are different from each other.
Two environments are considered equal if the distance of every neighboring atom differs by less than some tolerance.
The tolerance is specified in the *Tolerance (Å):* box.
In order to compare environments, all permutations of the atoms have to be considered and that calculation can be slow.

When the calculation has finished the app prints the number of environments that have been found and the average number of neighbors in the environments.

#### Units
* Cutoff should be given in angstrom
* Tolerance should be given in angstrom

### Analyze environments
The *Analyze environments* tab allows to visualize the calculated environments.
One can toggle between the *Unique* environments and *All* of them.
There is a slide bar to select the environment to visualize.

### Output environments
The *Output environments* tab prints the environments in Protein Data Bank (PDB) format.
One can toggle between the *Unique* environments and *All* of them.
The environments can be downloaded in a zip file.

## Acknowledgments

* The app uses several python libraries, for instance [ASE](https://wiki.fysik.dtu.dk/ase/), [NGLVIEW](https://github.com/arose/nglview), and [ipywidgets](https://ipywidgets.readthedocs.io/en/latest/index.html).
* I am grateful to Giovanni Pizzi and Dou Du for suggesting to deploy the tool using [Binder](https://mybinder.org/)+[appmode](https://github.com/oschuett/appmode).
* This tool was developed with support of the Swiss National Science Foundation (SNSF) through an Early Postdoc.Mobility fellowship.
* I also acknowledge funding from the NCCR MARVEL funded by the SNSF and from the CSI Computational Science Center funded by the Department of Energy of the USA.

## How to cite

If you are using this tool to find environments for enhanced sampling simulations please read and cite:
* [Pablo Piaggi, *Environment finder: A tool for finding and analyzing atomic environments in crystal structures*, Zenodo (2021)](https://doi.org/10.5281/zenodo.4746324)
* [Pablo Piaggi and Michele Parrinello, *Calculation of phase diagrams in the multithermal-multibaric ensemble*, J. Chem. Phys. 150, 244119 (2019)](https://aip.scitation.org/doi/full/10.1063/1.5102104)
