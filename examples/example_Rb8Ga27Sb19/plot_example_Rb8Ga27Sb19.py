from MCSPS.plotting import plot_swap_trajectory
from matplotlib import pyplot as plt

# Plot the annealed trajectories
plot_swap_trajectory('swaps*.out', min_energy=-0.22481528)
