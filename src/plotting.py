

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
  max_temp = -np.inf
  min_ene,max_ene = np.inf,-np.inf

  ach_min = []
  all_inds = []
  all_enes = []
  all_temps = []

  for fn in glob(file_pattern):

    inds,pinds,temps,enes = read_swap_trajectory(fn)
    mine,maxe = np.min(enes),np.max(enes)
    maxt = np.max(temps)
    maxi = np.max(inds)

    if maxi > max_ind:
      max_ind = maxi
    if mine < min_ene:
      min_ene = mine
    if maxe > max_ene:
      max_ene = maxe
    if maxt > max_temp:
      max_temp = maxt

    all_inds.append(inds)
    all_enes.append(enes)

    if min_energy is not None:
      t_mine = np.min(enes)
      if np.isclose(min_energy,t_mine,atol=1e-4) or t_mine < min_energy:
        ach_min.append(True)
      else:
        ach_min.append(False)

    if temperature:
      all_temps.append(temps)

  color_norm = plt.Normalize(0, max_temp)

  ordered_inds = list(range(len(all_inds)))
  if min_energy is not None:
    ordered_inds = [i for i,a in enumerate(ach_min) if not a]
    ordered_inds += [i for i,a in enumerate(ach_min) if a]

  for i in ordered_inds:

    inds = all_inds[i] + [max_ind]
    enes = all_enes[i] + [all_enes[i][-1]] 

    gcolors = colors[0]
    if min_energy is not None:
      if not ach_min[i]:
        gcolors = colors[1]

    # Color the segments by temperature
    if temperature:

      # Create the line segments
      points = np.array([inds,enes]).T.reshape(-1, 1, 2)
      segments = np.concatenate([points[:-1], points[1:]], axis=1)

      # Create the color mapping
      cmap = mpcol.LinearSegmentedColormap.from_list('', gcolors)

      lc = LineCollection(segments, cmap=cmap, norm=color_norm)
      lc.set_array(all_temps[i])
      col_line = ax.add_collection(lc)

    # Dont color the segments
    else:
      ax.plot(inds, enes, color=gcolors)

  # Set axis boundaries
  yoffset = 0.05 * (max_ene - min_ene)
  ax.set_xlim(0, max_ind*1.01)
  ax.set_ylim(min_ene-yoffset, max_ene+yoffset)

  if temperature:
    cbar = fig.colorbar(col_line)
    cbar.set_label('Temperature (K)')

  if show:
    plt.show()
    return None,None

  else:
    return fig,ax

