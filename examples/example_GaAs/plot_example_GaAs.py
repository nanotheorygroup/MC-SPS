from MCSPS.plotting import plot_swap_trajectory
from matplotlib import pyplot as plt

# Plot the absolute_zero trajectory
plot_swap_trajectory('swaps.out', temperature=False)

# Plot the annealed trajectory
plot_swap_trajectory('swaps_annealed.out')
