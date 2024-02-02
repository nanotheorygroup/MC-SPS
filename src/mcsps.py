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
  from .file_io import write_swap_accept
  from os.path import isfile
  from ase.io import write
  from ase import Atoms

  # Boltzmann const in eV
  kB = 1.68e-23/1.602e-19

  # Trajectory variables and structure information
  itr = 0
  last_swap_i = 0
  sswap = nswap
  nat = len(species)
  lattice = np.array(lattice)
  positions = np.array(positions)
  write_atoms = Atoms(species, positions=positions@lattice, cell=lattice, pbc=[1,1,1])

  # Verify that the input arrays have appropriate dimensions
  if len(lattice.shape) != 2 or not (lattice.shape[0] == 3 and lattice.shape[1] == 3):
    raise ValueError('Lattice shape must be (3,3)')
  if len(positions.shape) != 2 or positions.shape[0] != len(species) or positions.shape[1] != 3:
    raise ValueError('Positions shape must be (N,3), where N is the number of provided species')
  if nfixed > positions.shape[0]-2:
    raise ValueError(f'Cannot fix {nfixed} of {nat} sites. Decrease {nfixed} or provide more sites.')
  if len(set(species)) <= 1:
    raise ValueError('Species list must contain more than one type of species')

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
  write_swap_accept(swap_fname, 'w', itr, itr-last_swap_i, temperatures[0], ene)

  # Trajectory iteration
  for i,temp in enumerate(temperatures):
    for _ in range(temp_swaps[i]):
      itr += 1

      rswap = 1 + np.random.randint(np.min((nat-nfixed,nswap)))
      t_species = sps_swap(species, nfixed=nfixed, nswaps=rswap)
      t_ene = calculator.predict_formation_energy(lattice, t_species, positions)

      dE = t_ene - ene
      boltz = False if temp==0 else np.exp(-dE/(kB*temp)) > np.random.rand()

      # Accept condition
      if dE < 0 or boltz:
        write_swap_accept(swap_fname, 'a', itr, itr-last_swap_i, temp, t_ene)
        ene = t_ene
        nswap = sswap
        last_swap_i = itr
        species = t_species
        if ene < emin:
          emin = ene
          if emin_fname is not None:
            write_atoms.set_chemical_symbols(species)
            write(emin_fname, write_atoms)

      else:
        if (itr-last_swap_i) % nswap_inc == 0:
          nswap += 1

      if stop is not None:
        if ene <= stop:
          return



def sps_vacancy ( lattice,
                  vacancy_species,
                  vacancy_positions,
                  temperatures,
                  temp_swaps=[],
                  fixed_species=None,
                  fixed_positions=None,
                  nswap = 2,
                  nswap_inc = 1000,
                  swap_fname = 'swaps.out',
                  emin_fname = 'structure.emin.xyz',
                  occupation_factors=None,
                  calculator = None,
                  stop=None ):
  '''
    Perform the SPS routine on a set of atomic sites that contain some vacant sites.

    Arguments:
      lattice (list or ndarray): 3x3 matrix representing the three lattice vectors [R1, R2, R3]
      fixed_species (list): List of atomic symbols for sites that cannot participate in swapping
      fixed_positions (list or ndarray): Nx3 matrix, with the crystal 3-coordinate for each of N fixed atomic sites
      vacancy_species (list): List of atomic symbols for sites participating in the swaping
      vacancy_sites (list or ndarray): Nx3 matrix, with crystal 3-coordinate of sites that can be occupied or unoccupied by vacancy_species
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
  from .file_io import write_swap_accept
  from os.path import isfile
  from ase.io import write
  from ase import Atoms

  # Boltzmann const in eV
  kB = 1.68e-23/1.602e-19

  # Trajectory variables
  itr = 0
  last_swap_i = 0
  sswap = nswap

  # Set up fixed positions
  if None in [fixed_species, fixed_positions]:
    if fixed_species is not None or fixed_positions is not None:
      print('User must provide 1:1::species:position to fix. No fixed species or positions will be included.')
    fixed_species = []
    fixed_positions = np.empty((0,3), dtype=float)
  else:
    fixed_positions = np.array(fixed_positions)

  # Set up lattice, species, vacancy positions
  species = fixed_species + vacancy_species
  lattice = np.array(lattice)
  vacancy_positions = np.array(vacancy_positions)

  # Verify that the input arrays have appropriate dimensions
  if len(lattice.shape) != 2 or not (lattice.shape[0] == 3 and lattice.shape[1] == 3):
    raise ValueError('Lattice shape must be (3,3)')
  if len(fixed_positions.shape) != 2 or fixed_positions.shape[0] != len(fixed_species) or fixed_positions.shape[1] != 3:
    raise ValueError('Fixed positions shape must be (N,3), where N is the number of fixed species')
  if len(vacancy_positions.shape) != 2 or vacancy_positions.shape[0] < len(vacancy_species):
    msg = f'Vacancy positions shape must be (N,3), where N is greater than or equal to the number of vacancy species.'
    raise ValueError(msg)

  # Static array lengths
  nvspecies = len(vacancy_species)
  nfixed = fixed_positions.shape[0]
  nvsites = vacancy_positions.shape[0]

  # Set the initial vacancy occupation
  vacancy_indices = [i for i in range(nvspecies)]
  vacancy_occupation = {i:i for i in vacancy_indices}
  v_positions = [vacancy_positions[i] for i in vacancy_occupation]

  # Create Atoms object for output
  positions = np.concatenate([fixed_positions, v_positions])
  write_atoms = Atoms(species, positions=positions@lattice, cell=lattice, pbc=[1,1,1])

  # If there is no provided calculator, initialize the default M3GNet calculator
  if calculator is None:
    from .calculators import MEGNet_Calculator
    calculator = MEGNet_Calculator()

  # If the number of swaps per temperature is undefined, assign each to 1
  if len(temp_swaps) == 0:
    temp_swaps = np.ones(len(temperatures), dtype=np.short)
  if len(temperatures) != len(temp_swaps):
    raise ValueError('temperatures and temp_swaps must contain the same number of elements.')

  ene = emin = calculator.predict_formation_energy(lattice, species, positions)

  # Initialize the swaps output file. Exit if the file exists already.
  if isfile(swap_fname):
    raise FileExistsError(f'File {swap_fname} already exists. Will not overwrite.')
  write_swap_accept(swap_fname, 'w', itr, itr-last_swap_i, temperatures[0], ene, occupation_factors=occupation_factors)

  # Trajectory iteration
  for i,temp in enumerate(temperatures):
    for _ in range(temp_swaps[i]):
      itr += 1

      tvinds = vacancy_indices.copy()
      tvocc = vacancy_occupation.copy()

      # Perform stochastic swaps
      ns, nc = 0, 1+np.random.randint(nswap)
      ri = np.random.randint(nvspecies)
      while ns < nc:
        nri = np.random.randint(nvsites)
        if nri in tvocc:
          if vacancy_species[tvocc[nri]] == vacancy_species[ri]:
            continue
          tri = tvocc[nri]
          tvocc[tvinds[ri]],tvocc[nri] = tvocc[nri],tvocc[tvinds[ri]]
          tvinds[ri],tvinds[tri] = nri,tvinds[ri]
        else:
          del tvocc[tvinds[ri]]
          tvinds[ri] = nri
          tvocc[nri] = ri
        ns += 1

      # Update positions arrays
      v_positions = np.array([vacancy_positions[i] for i in tvinds])
      positions[nfixed:,:] = v_positions

      # Calculate energy and evaluate Metropolis condition
      t_ene = calculator.predict_formation_energy(lattice, species, positions)
      dE = t_ene - ene
      boltz = False if temp==0 else np.exp(-dE/(kB*temp)) > np.random.rand()
      if dE < 0 or boltz:
        write_swap_accept(swap_fname, 'a', itr, itr-last_swap_i, temp, t_ene, occupation_factors=occupation_factors)
        ene = t_ene
        nswap = sswap
        last_swap_i = itr
        vacancy_indices = tvinds
        vacancy_occupation = tvocc
        if ene < emin:
          emin = ene
          if emin_fname is not None:
            write_atoms.set_scaled_positions(positions)
            write(emin_fname, write_atoms)

      else:
        # Increase the number of swaps after nswap_inc rejections
        if (itr-last_swap_i) % nswap_inc == 0:
          nswap += 1

      # Halt if stop condition is met
      if stop is not None:
        if ene <= stop:
          return
