
# Monte Carlo Site Permutation Search

Toolkit for optimizing atomic decoration on fixed geometries. MEGNet or M3GNet provide rapid estimation of total energy, which can accurately rank the energetic favorability of various configurations. Simulated annealing and basin hopping techniques assist in discovery of locally and globally minimal configurations.  
Functionality and results of this software is documented at https://doi.org/10.21203/rs.3.rs-3740642/v1. The work is currently under review with NPJ Computational Materials.  
  
## Requirements:
  * ASE
  * Numpy
  * MEGNet or M3GNet  


## Installation:
  From inside of the MC-SPS root directory, run the following python command.
  * python setup.py install


## Usage:
  * The following SPS routines can be imported from the mcsps module, sps\_fixed, sps\_vacancy, sps\_cluster. Structure lattice, atomic basis, and temperature trajectory are supplied directly to the SPS routines.
  * Examples documenting the package usage are located in the examples directory.
  * The main\*.py scripts can be run in the background, and output can be monitored through the trajectory output file and structure output file. Trajectory output is updated upon each accepted configuration. Structure output is updated each time a structure is identified with a lower total energy.
  * Visualization tools and advanced features for documenting site occupation factors and saving multiple favorable structures are coming soon. 


## Examples:
  * GaAs - 3x3x3 Gallium Arsenide supercell. Demonstrates the most basic features of the SPS package.
  * Rb8Ga27Sb19 - Superstructural ordering predicted on the elemental type-I clathrate geometry
  * Vacancy ordered SPS optimization examples coming soon. 

