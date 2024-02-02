from MCSPS.utilities import create_supercell
from MCSPS.mcsps import sps_fixed
import numpy as np

# Seed the trajectory for reproducibility
np.random.seed(4321)

# Create the GaAs unit lattice and atomic basis
fcc = 0.5 * np.array([[-1,0,1],[0,1,1],[-1,1,0]], dtype=float)
unit_lattice = 5.75 * fcc

unit_species = ['Ga', 'As']
unit_positions = np.array([[0,0,0], [0.25,0.25,0.25]])

# Create a 3x3x3 supercell
supercell_dimensions = [4, 4, 4]
lattice, positions, species = create_supercell(unit_lattice, unit_positions, unit_species, supercell_dimensions)

# Randomize the list of species
np.random.shuffle(species)

# Define temperature trajectory and number of steps per temperature
#  This trajectory runs from 100K to 0K in 200 steps, performing 5 swaps at each step
temperatures = np.linspace(80, 0, 201)
num_swaps_per_temperature = 200 * [20] + [2000]

# Start SPS. The output filenames are specified explicitly
sps_fixed(lattice, species, positions, temperatures, num_swaps_per_temperature,
          swap_fname='swaps_annealed.out', emin_fname='structure_annealed.emin.xyz')

