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
* Two algorithms for finding unique environments are available. A sorting algorithm that is fast but in some cases can be wrong. A permutation algorithm that is slow but always yield correct results. My advice is to use the fast version and change to the slow one if the algorithm finds more unique environments than expected.
* Visualize the environments.
* Output the environments in pdb format.

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


## Acknowledgments

The app uses several python libraries, for instance ASE and NGLVIEW.
