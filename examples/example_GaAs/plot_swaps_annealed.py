from matplotlib.collections import LineCollection
from matplotlib import colors as mpcol
from matplotlib import pyplot as plt
import numpy as np

from MCSPS.file_io import read_swap_trajectory

swap_fname = 'swaps_annealed.out'

inds,pinds,temps,enes = read_swap_trajectory(swap_fname)

fig,ax = plt.subplots(figsize=(6,4))
fig.suptitle('GaAs Swap Trajectory')
ax.set_xlabel('Iteration')
ax.set_ylabel('Energy (eV/atom)')

# Create the line segments
points = np.array([inds,enes]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Create the color mapping
norm = plt.Normalize(0, np.amax(temps))
cmap = mpcol.LinearSegmentedColormap.from_list('', ['blue','red'])

lc = LineCollection(segments, cmap=cmap, norm=norm)
lc.set_array(temps)
line = ax.add_collection(lc)

ax.set_xlim(0,1000)
ax.set_ylim(-1,1)
#ax.plot(inds, enes, color='black')

plt.show()

