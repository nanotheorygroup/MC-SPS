


def create_supercell ( lattice, species, positions, nsc=(2,2,2) ):
  '''
    Create a supercell lattice and set of atomic positions from a unit arrangement
  '''
  import numpy as np

  supercell_species = []
  supercell_positions = []

  for i in range(nsc[0]):
    for j in range(nsc[1]):
      for k in range(nsc[2]):
        origin = np.array([i,j,k])
        for spec,pos in zip(species, positions):
          supercell_species.append(spec)
          supercell_positions.append(pos + origin)

  supercell_lattice = lattice.copy()
  supercell_positions = np.array(supercell_positions)
  for i in range(3):
    supercell_lattice[i] *= nsc[i]
    supercell_positions[:,i] /= nsc[i]
     
  return supercell_lattice,supercell_species,supercell_positions



def boundary_verification ( lattice, species, positions, nfixed, nvacant=0 ):

  if len(lattice.shape) != 2 or not (lattice.shape[0] == 3 and lattice.shape[1] == 3):
    raise ValueError('Lattice shape must be (3,3)')

  if len(lattice.shape) != 2 or positions.shape[0] != len(species) or positions.shape[1] != 3:
    raise ValueError('Positions shape must be (N,3), where N is the number of provided Species')

  if nfixed > positions.shape[0]-2:
    raise ValueError(f'Cannot fix {nfixed} of {nat} sites. Decrease {nfixed} or provide more sites.')

  if len(set(species)) <= 1:
    raise ValueError('Species list must contain more than one type of species')

# Check nfixed + vacancies == nspecies

