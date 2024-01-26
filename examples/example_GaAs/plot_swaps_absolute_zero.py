from MCSPS.file_io import read_swap_trajectory
from matplotlib import pyplot as plt

swap_fname = 'swaps.out'

inds,pinds,temps,enes = read_swap_trajectory(swap_fname)

fig,ax = plt.subplots(figsize=(6,4))
fig.suptitle('GaAs Swap Trajectory')
ax.set_xlabel('Iteration')
ax.set_ylabel('Energy (eV/atom)')

ax.plot(inds, enes, color='black')

plt.show()

