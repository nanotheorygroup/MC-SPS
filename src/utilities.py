


def create_supercell ( lattice, positions, species=None, nsc=(2,2,2) ):
  '''
    Create a supercell lattice and set of atomic positions from a unit arrangement.
    If species are provided, they will be duplicated accordingly.

    Arguments:
      lattice (ndarray): 3x3 matrix with the lattice vectors
      positions (ndarray): Nx3 matrix containing the unit positions as 3-vectors
      species (list): List of species to duplicate. If None, only lattice and positions are replicated
      nsc (tuple, list, or ndarray): Number of cells to replicate in each direction (nx,ny,nz)

    Returns:
      (ndarray,ndarray): supercell lattice and positions, if species is not specified
      (ndarray,ndarray,list): supercell lattice, positions, and replicated species, if species is specified
  '''
  import numpy as np

  supercell_positions = []

  for i in range(nsc[0]):
    for j in range(nsc[1]):
      for k in range(nsc[2]):
        origin = np.array([i,j,k])
        supercell_positions += [pos+origin for pos in positions]

  supercell_lattice = lattice.copy()
  supercell_positions = np.array(supercell_positions)
  for i in range(3):
    supercell_lattice[i] *= nsc[i]
    supercell_positions[:,i] /= nsc[i]

  if species is not None:
    supercell_species = np.prod(nsc) * species
    return supercell_lattice,supercell_positions,supercell_species

  else:
    return supercell_lattice,supercell_positions

