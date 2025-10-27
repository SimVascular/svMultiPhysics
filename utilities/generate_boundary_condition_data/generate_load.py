import numpy as np
import os

# Go to directory of this script
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Number of timesteps and number of Fourier modes
n_timesteps = 101
n_modes = 64

# Generate time values from 0 to 2
time = np.linspace(0, 2, n_timesteps)

# Generate ramp from 0 to 10 in 0.5 seconds, then hold until 2 seconds.
load = np.zeros(n_timesteps)
load[time < 0.5] = 10 * time[time < 0.5] / 0.5
load[(time >= 0.5)] = 10


# Write the time and stress values to a text file
with open("load.dat", "w") as file:
    file.write(f"{n_timesteps} {n_modes}\n")
    for t, s in zip(time, load):
        file.write(f"{t:.3f} {s:.3f}\n")

# Plot the stress values
import matplotlib.pyplot as plt
plt.plot(time, load)
plt.xlabel("Time (s)")
plt.ylabel("Load (dynes/cm^2)")
plt.title("Load")
plt.savefig("load.png")