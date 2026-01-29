import numpy as np
import os

# Go to directory of this script
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Number of timesteps and number of Fourier modes
n_timesteps = 1001
n_modes = 256

# Generate time values from 0 to 1.5. Go to 1.5 to avoid Gibbs phenomenon.
time = np.linspace(0, 1.5, n_timesteps)

# Generate ramp stress values with a maximum of 50
max_stress = 50
stress = max_stress * time  # Simple linear ramp from 0 to max_stress

# Write the time and stress values to a text file
with open("stress.dat", "w") as file:
    file.write(f"{n_timesteps} {n_modes}\n")
    for t, s in zip(time, stress):
        file.write(f"{t:.3f} {s:.3f}\n")

# Plot the stress values
import matplotlib.pyplot as plt
plt.plot(time, stress)
plt.xlabel("Time (s)")
plt.ylabel("Stress (dynes/cm^2)")
plt.title("Active stress (ramp)")
plt.savefig("stress.png")