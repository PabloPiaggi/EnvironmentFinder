# Environment finder

A tool for finding and analyzing atomic environments in crystal structures.

Launch the app! [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PabloPiaggi/EnvironmentFinder/master?urlpath=apps%2FEnvironmentFinder.ipynb)

## Purpose

The main purpose of this app is to find atomic environments in crystal structures.
These environments can then be used to define measures of similarity and collective variables for enhanced sampling simulations.
The app can also be used to visualize chemical environments around an atom for general purposes.

## Features

* Choose examples or upload own configuration.
* Determine unique environments in crystal structures. Generally, these environments correspond to the *basis* of the crystal structure.
* Two algorithms for finding unique environments are available. 
	* A sorting algorithm that is fast but in some pathological cases can be wrong. 
	* A permutation algorithm that is slow but always yield correct results. My advice is to use the fast version and change to the slow one if the algorithm finds more unique environments than expected.
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
The *Choose configuration* tab has two toggle buttons *Choose* and *Upload*.
By default the *Choose* option allows the user to select from some examples.
These examples are useful to play with and get a grasp of the tool.
The *Upload* option should be used to analyze a configuration stored in the user's hard drive.
Once a configuration file has been uploaded the user can switch to the *Choose* option and the file should appear in the dropdown choices.
Any configuration supported by the ASE library can be used.
A complete list can be found [here](https://wiki.fysik.dtu.dk/ase/ase/io/io.html).
ASE determines the file type by its extension so be sure to name the file properly!

### Define environments
The *Define environments* tab is where the user defines the environments.
There are three different ways to define the environments, *Type*, *String*, and *Step* that can be selected with the toggle buttons.
Typically the most useful one is the *Type* option.
In this case one chooses the atom type of the central atom (A), the atom type of the neighboring atoms (B), and a cutoff.
These parameters are set using the *Atom type 1:*, *Atom type 2:*, and *Cutoff:* boxes.
The app will find all environments with atoms of type B around atoms of type A inside the cutoff.

There are also two checkboxes *Find unique environments?* and *Fast algorithm*.
The *Find unique environments?* checkbox triggers the calculation of unique environments.
For this purpose it compares all environments and keeps only the ones that are different from each other.
Two environments are considered equal if the distance of every neighboring atom differs by less than some tolerance.
The tolerance is specified in the *Tolerance:* box.
In order to compare environments, all permutations of the atoms have to be considered and that calculation can be slow.
The *Fast algorithm* checkbox switches between a sorting algorithm that is fast but in some pathological cases can be wrong, and a permutation algorithm that is slow but always yield correct results.
My advice is to use the fast version and change to the slow one if the algorithm finds more unique environments than expected.

When the calculation has finished the app prints the number of environments that have been found.

#### Units
* Cutoff should be given in angstrom
* Tolerance should be given in angstrom

### Analyze environments
The *Analyze environments* tab allows to visualize the calculated environments.
One can toggle between the *Unique* environments and *All* of them. 
Bear in mind that you might need to toggle between these options to refresh the results.

There is a slide bar to select the environment to visualize.

### Output environments
The *Output environments* tab prints the environments in Protein Data Bank (PDB) format.
One can toggle between the *Unique* environments and *All* of them. 
Bear in mind that you might need to toggle between these options to refresh the results.

There is a slide bar to select the environment to visualize.

## Acknowledgments

The app uses several python libraries, for instance [ASE](https://wiki.fysik.dtu.dk/ase/) and [NGLVIEW](https://github.com/arose/nglview).
I am grateful to Giovanni Pizzi and Dou Du for suggesting to deploy the tool using [Binder](https://mybinder.org/)+[appmode](https://github.com/oschuett/appmode).
This tool was developed with support of the Swiss National Science Foundation (SNSF) through an Early Postdoc.Mobility fellowship.
I also acknowledge funding from the NCCR MARVEL funded by the SNSF and from the CSI Computational Science Center funded by the Department of Energy of the USA.
