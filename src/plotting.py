

def plot_swap_trajectory ( file_pattern, temperature=True, min_energy=None, show=True ):
  from matplotlib.collections import LineCollection
  from matplotlib import colors as mpcol
  from matplotlib import pyplot as plt
  from glob import glob
  import numpy as np

  from .file_io import read_swap_trajectory

  fig,ax = plt.subplots(figsize=(6,4))
  fig.suptitle('GaAs Swap Trajectory')
  ax.set_xlabel('Iteration')
  ax.set_ylabel('Energy (eV/atom)')


  colors = ['black', 'gray']
  col_line = None
  if temperature:
    colors = [['blue','red'], ['lightblue','pink']]

  max_ind = -np.inf
  min_ene,max_ene = np.inf,-np.inf

  for fn in glob(file_pattern):

    inds,pinds,temps,enes = read_swap_trajectory(fn)
    mine,maxe = np.min(enes),np.max(enes)
    maxi = np.max(inds)

    if maxi > max_ind:
      max_ind = maxi
    if mine < min_ene:
      min_ene = mine
    if maxe > max_ene:
      max_ene = maxe

    gcolors = colors[0]
    if min_energy is not None and not np.isclose(min_energy,np.min(enes)):
      gcolors = colors[1]

    # Color the segments by temperature
    if temperature:

      # Create the line segments
      points = np.array([inds,enes]).T.reshape(-1, 1, 2)
      segments = np.concatenate([points[:-1], points[1:]], axis=1)

      # Create the color mapping
      norm = plt.Normalize(0, np.amax(temps))
      cmap = mpcol.LinearSegmentedColormap.from_list('', gcolors)

      lc = LineCollection(segments, cmap=cmap, norm=norm)
      lc.set_array(temps)
      col_line = ax.add_collection(lc)

    # Dont color the segments
    else:
      ax.plot(inds, enes, color=gcolors)

  yoffset = 0.05 * (max_ene - min_ene)
  ax.set_xlim(0, max_ind*1.1)
  ax.set_ylim(min_ene-yoffset, max_ene+yoffset)

  if temperature:
    cbar = fig.colorbar(col_line)
    cbar.set_label('Temperature (K)')

  if show:
    plt.show()
    return None,None

  else:
    return fig,ax

