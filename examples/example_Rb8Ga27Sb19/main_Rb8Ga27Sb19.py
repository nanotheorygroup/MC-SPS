from MCSPS.mcsps import sps_fixed
from ase.io import read
import numpy as np

# Seed the trajectory for reproducibility
np.random.seed(1234)

# Read the clathrate structure from file
atoms = read('Rb32Si184.xyz')

lattice = atoms.cell[:]
species = atoms.get_chemical_symbols()
positions = atoms.positions

# Define the number of fixed atoms, at the front of the species list.
#  Here, the Rubidium atoms will not undergo exchange
n_fixed_atoms = 32

# Set the remaining species list to the correct stoichiometry
species[n_fixed_atoms:] = 4 * (27*['Ga'] + 19*['Sb'])

# Define temperature trajectory and number of steps per temperature
#  This trajectory runs from 30K to 0K in 301 steps, performing 1000 swaps at each temperature
temperatures = np.linspace(30, 0, 301)
num_swaps_per_temperature = 301 * [1000]

# Run 50 trials
for i in range(50):

  swap_fname = f'swaps.{i}.out'
  emin_fname = f'structure.{i}.xyz'

  # Start SPS
  sps_fixed(lattice, species, positions, temperatures,
            temp_swaps = num_swaps_per_temperature,
            nfixed = n_fixed_atoms,
            swap_fname = swap_fname,
            emin_fname = emin_fname)

