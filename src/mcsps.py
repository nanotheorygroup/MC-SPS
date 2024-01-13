import numpy as np


def sps_swap ( species:list, nfixed=0, nswaps=1 ) -> dict:
  '''
    Interchange pairs in a list, stochastically

    Arguments:
      species (list): The list of elements to interchange
      nfixed (int): The number of elements to remain fixed, which must be placed at the beginning of the list
      nswaps (int): The number of swaps to perform

    Returns:
      (list): The resulting list with elements exchanged
  '''

  nat = len(species)
  t_species = species.copy()

  for _ in range(nswaps):

    ri1 = np.random.randint(nfixed, nat)
    ri2 = np.random.randint(nfixed, nat)

    while ri1 == ri2 or species[ri1] == species[ri2]:
      ri2 = np.random.randint(nfixed, nat)

    t_species[ri1],t_species[ri2] = t_species[ri2],t_species[ri1]

  return t_species



def sps_fixed ( lattice,
                species,
                positions,
                temperatures,
                temp_swaps=[],
                nfixed = 0,
                nswap = 2,
                nswap_inc = 1000,
                swap_fname = 'swaps.out',
                emin_fname = 'structure.emin.xyz',
                calculator = None,
                stop=None ):
  '''
    Perform the SPS routine on a fixed atomic basis without vacant sites.

    Arguments:
      lattice (list or ndarray): 3x3 matrix representing the three lattice vectors [R1, R2, R3]
      species (list): List of atomic symbols for each constituent site
      positions (list or ndarray): Nx3 matrix, with the crystal 3-coordinate for each of N atomic sites
      temperatures (list or ndarray): Temperature trajectory for simulated annealing trials
      temp_swaps (list or ndarray): Number of swaps to perform at each temperature in the trajectory. (Must contain one number for each provided temperature.
      nfixed (int): Number of sites to neglect from the swapping routine. Fixed sites must come first in the species and positions lists.
      nswap (int): Number of swaps to be performed at each step. This value is increased after nswap_inc rejected iterations.
      nswap_inc (int): Number of rejected trial configurations performed before increasing nswap.
      swap_fname (str): File name for the trajectory output file.
      emin_filename (str): File name for the minimum energy structure, output in the xyz format.
      calculator (str): Calculator for evaluating the structure energy.
      stop (float): Terminate the trial if the predicted structure energy is at or below the provided stop value.
  '''
  from .utilities import boundary_verification
  from os.path import isfile
  from ase.io import write
  from ase import Atoms

  kB = 1.68e-23/1.602e-19 # Boltzmann const in eV

  itr = 0
  last_swap_i = 0
  nswap = sswap = 2
  nat = len(species)
  lattice = np.array(lattice)
  positions = np.array(positions)
  write_atoms = Atoms(species, positions=positions, cell=lattice, pbc=[1,1,1])

  # Verify that the input arrays have appropriate dimensions
  boundary_verification(lattice, species, positions, nfixed)

  # If the number of swaps per temperature is undefined, assign each to 1
  if len(temp_swaps) == 0:
    temp_swaps = np.ones(len(temperatures), dtype=np.short)

  if len(temperatures) != len(temp_swaps):
    raise ValueError('temperatures and temp_swaps must contain the same number of elements.')

  # If there is no provided calculator, initialize the default M3GNet calculator
  if calculator is None:
    from .calculators import MEGNet_Calculator
    calculator = MEGNet_Calculator()

  ene = emin = calculator.predict_formation_energy(lattice, species, positions)

  # Initialize the swaps output file. Exit if the file exists already.
  if isfile(swap_fname):
    raise FileExistsError(f'File {swap_fname} already exists. Will not overwrite.')
  with open(swap_fname, 'w') as f:
    f.write('{} {} {} {}\n'.format(itr, itr-last_swap_i, temperatures[0], ene))

  for i,temp in enumerate(temperatures):
    for _ in range(temp_swaps[i]):
      itr += 1

      rswap = 1 + np.random.randint(np.min((nat-nfixed,nswap)))
      t_species = sps_swap(species, nfixed=nfixed, nswaps=rswap)
      t_ene = calculator.predict_formation_energy(lattice, t_species, positions)

      dE = t_ene - ene
      boltz = False if temp==0 else np.exp(-dE/(kB*temp)) > np.random.rand()
      if dE < 0 or boltz:
        with open(swap_fname, 'a') as f:
          f.write('{} {} {} {}\n'.format(itr, itr-last_swap_i, temp, t_ene))
        ene = t_ene
        nswap = sswap
        last_swap_i = itr
        species = t_species
        if ene < emin:
          emin = ene
          write_atoms.set_chemical_symbols(species)
          write(emin_fname, write_atoms)
          # Write XYZ structure file for minimum energy structure

      else:
        if (itr-last_swap_i) % nswap_inc == 0:
          nswap += 1

      if stop is not None:
        if ene <= stop:
          return
