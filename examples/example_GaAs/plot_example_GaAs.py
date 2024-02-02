from MCSPS.plotting import plot_swap_trajectory

# Plot the absolute_zero trajectory
plot_swap_trajectory('swaps.out', temperature=False, show=False)

# Plot the annealed trajectory
plot_swap_trajectory('swaps_annealed.out')
