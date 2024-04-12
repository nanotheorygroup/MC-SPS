from MCSPS.utilities import create_supercell
from multiprocessing import Process,Manager
from MCSPS.mcsps import sps_fixed
import numpy as np

# Probe the order-disorder phase transition in CuZn
#  Trajectories are performed in parallel

# Create the GaAs unit lattice and atomic basis
unit_lattice = 2.955 * np.eye(3)

unit_species = ['Cu', 'Zn']
unit_positions = np.array([[0,0,0], [0.5,0.5,0.5]])

# Create a 5x5x5 supercell
supercell_dimensions = [5, 5, 5]
lattice, positions, species = create_supercell(unit_lattice, unit_positions, unit_species, supercell_dimensions)

# Temperature trajectory
ntemp = 12
nstep = 25000
nswaps = ntemp * [nstep]
temps = np.linspace(4.0, 1.8, ntemp)

# Start SPS with multiprocessing
def start_sps ( fn_swap, fn_smin ):
  sps_fixed(lattice, species, positions, temps, nswaps,
            swap_fname=fn_swap, emin_fname=fn_smin)

# Start multiprocess manager
manager = Manager()
dos = manager.dict()

# Set number of processes to start
nproc = 4
procs = []

# Create the annealing processes
for i in range(nproc):

  np.random.seed(i*4321)
  fname = f'swaps.{i}.out'
  fn_smin = f'structure.emin.{i}.txt'

  p = Process(target=start_sps, args=(fname, fn_smin))
  procs.append(p)
  procs[-1].start()

# Join processes
for p in procs:
  p.join()
