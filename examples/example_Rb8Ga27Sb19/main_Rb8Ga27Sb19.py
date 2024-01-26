from MCSPS.mcsps import sps_fixed
from ase.io import read
import numpy as np

# Seed the trajectory for reproducibility
np.random.seed(54321)

# Read the clathrate structure from file
atoms = read('Rb32Si184.xyz')

lattice = atoms.cell[:]
species = atoms.get_chemical_symbols()
positions = atoms.get_scaled_positions()

# Define the number of fixed atoms, at the front of the species list.
#  Here, the Rubidium atoms will not undergo exchange
n_fixed_atoms = 32

# Set the remaining species list to the correct stoichiometry on arbitrary sites
species[n_fixed_atoms:] = 4 * (27*['Ga'] + 19*['Sb'])

# Define temperature trajectory and number of steps per temperature
#  This trajectory runs from 30K to 0K in 301 steps, performing 1000 swaps at each temperature
temperatures = np.linspace(30, 0, 211)
num_swaps_per_temperature = 210 * [1000] + [10000]

# Run 5 trials
for i in range(5):

  swap_fname = f'swaps.{i}.out'
  emin_fname = f'structure.{i}.xyz'

  # Start SPS
  sps_fixed(lattice, species, positions, temperatures,
            temp_swaps = num_swaps_per_temperature,
            nfixed = n_fixed_atoms,
            nswap_inc = 250,
            swap_fname = swap_fname,
            emin_fname = emin_fname)

