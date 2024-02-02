from MCSPS.utilities import create_supercell
from MCSPS.mcsps import sps_vacancy
from ase.io import read
import numpy as np

np.random.seed(12345)

# Perovskite unit lattice and atomic basis
cubic = 6.29 * np.eye(3, dtype=float)
unit_positions = 0.5 * np.array([[0,0,0],[1,1,1],[1,1,0],[1,0,1],[0,1,1]], dtype=float)

# Create supercell
supercell_dimensions = (2, 2, 2)
lattice, positions = create_supercell(cubic, unit_positions, nsc=supercell_dimensions)

# Assign species, omitting 4 Sn for the vacant sites
species = 4 * (2*['Cs'] + ['Sn'] + 6*['I'])

# Perform 3000 swaps at 0K
trange = [0]
srange = [3000]

# Perform the SPS optimization 10 times
for i in range(10):
  swap_fname = f'swaps_0K.{i}.out'
  structure_fname = f'structure_0K.{i}.xyz'
  sps_vacancy(lattice, species, positions, trange, srange, 
              swap_fname=swap_fname, emin_fname=structure_fname)

