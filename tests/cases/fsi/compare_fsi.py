"""
Compare monolithic vs partitioned FSI results and generate:
1. Field comparison at final time step
2. Videos of both simulations
3. Coupling performance plot for partitioned FSI
"""

import numpy as np
import meshio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.animation as animation
import os, re

mono_dir = "pipe_3d/1-procs"
part_dir = "pipe_3d_partitioned/1-procs"
coupling_log = "pipe_3d_partitioned/coupling_log.txt"
out_dir = "pipe_3d_partitioned"

n_steps = 50

# ============================================================
# 1. Parse coupling log
# ============================================================
print("Parsing coupling log...")
coupling_data = {}  # ts -> [(outer, omega, rel_disp), ...]

with open(coupling_log) as f:
    for line in f:
        # Format: CP TS-ITER  TIME  [dB REL_DISP OMEGA]
        m = re.match(r'\s*CP\s+(\d+)-(\d+)\s+([\d.e+-]+)\s+\[(-?\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)\]', line)
        if m:
            ts = int(m.group(1))
            outer = int(m.group(2))
            rel_disp = float(m.group(5))
            omega = float(m.group(6))
            if ts not in coupling_data:
                coupling_data[ts] = []
            coupling_data[ts].append((outer, omega, rel_disp))

# ============================================================
# 2. Coupling performance plot
# ============================================================
print("Creating coupling performance plot...")

fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# Top: number of coupling iterations per time step
ts_list = sorted(coupling_data.keys())
n_iters = [len(coupling_data[ts]) for ts in ts_list]
axes[0].bar(ts_list, n_iters, color='steelblue', width=0.8)
axes[0].set_ylabel('Coupling iterations')
axes[0].set_title('Partitioned FSI Coupling Performance (Aitken relaxation)')
axes[0].set_ylim(0, max(n_iters) + 1)

# Bottom: convergence history (all time steps overlaid)
cmap = plt.cm.viridis
for i, ts in enumerate(ts_list):
    iters = coupling_data[ts]
    x = [it[0] for it in iters]
    y = [it[2] for it in iters]
    color = cmap(i / max(len(ts_list) - 1, 1))
    axes[1].semilogy(x, y, 'o-', color=color, alpha=0.6, lw=1, markersize=4)

# Highlight first and last
for ts, color, label in [(ts_list[0], 'blue', f'TS {ts_list[0]}'),
                          (ts_list[-1], 'red', f'TS {ts_list[-1]}')]:
    iters = coupling_data[ts]
    x = [it[0] for it in iters]
    y = [it[2] for it in iters]
    axes[1].semilogy(x, y, 'o-', color=color, lw=2, markersize=6, label=label)

axes[1].axhline(1e-8, color='green', linestyle='--', linewidth=1.5, label='Tolerance (1e-8)')
axes[1].set_xlabel('Coupling iteration within time step')
axes[1].set_ylabel('Relative displacement change')
axes[1].legend(loc='upper right')
axes[1].set_ylim(1e-12, 10)
axes[1].set_xlim(0.5, max(n_iters) + 0.5)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'coupling_performance.png'), dpi=150)
plt.close()
print(f"  Saved {out_dir}/coupling_performance.png")

# ============================================================
# 3. Final time step comparison
# ============================================================
print("Comparing final time step fields...")

mono = meshio.read(os.path.join(mono_dir, f"result_{n_steps:03d}.vtu"))
part = meshio.read(os.path.join(part_dir, f"result_{n_steps:03d}.vtu"))

print(f"\n{'Field':<20} {'Mono max':>12} {'Part max':>12} {'Rel diff':>12}")
print("-" * 60)
for f in mono.point_data:
    m = np.max(np.abs(mono.point_data[f]))
    if f in part.point_data:
        p = np.max(np.abs(part.point_data[f]))
    else:
        p = 0.0
    rd = abs(m - p) / (m + 1e-30)
    print(f"{f:<20} {m:12.4e} {p:12.4e} {rd:12.4e}")

# ============================================================
# 4. Videos (side-by-side velocity magnitude on z=0 slice)
# ============================================================
print("\nCreating videos...")

def get_velocity_mag(fname):
    """Read VTU and return (points, velocity_magnitude)"""
    m = meshio.read(fname)
    pts = m.points
    vel = m.point_data.get('Velocity', np.zeros((len(pts), 3)))
    vmag = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2 + vel[:, 2]**2)
    return pts, vmag

def get_displacement_mag(fname):
    """Read VTU and return (points, displacement_magnitude)"""
    m = meshio.read(fname)
    pts = m.points
    # Try FS_Displacement first (monolithic), then Displacement
    for key in ['FS_Displacement', 'Displacement']:
        if key in m.point_data:
            d = m.point_data[key]
            if np.max(np.abs(d)) > 1e-15:
                return pts, np.sqrt(d[:, 0]**2 + d[:, 1]**2 + d[:, 2]**2)
    return pts, np.zeros(len(pts))

# Get all data for velocity video
print("  Reading velocity data...")
mono_vel = []
part_vel = []
for ts in range(1, n_steps + 1):
    mono_pts, mono_vmag = get_velocity_mag(os.path.join(mono_dir, f"result_{ts:03d}.vtu"))
    part_pts, part_vmag = get_velocity_mag(os.path.join(part_dir, f"result_{ts:03d}.vtu"))
    mono_vel.append((mono_pts, mono_vmag))
    part_vel.append((part_pts, part_vmag))

# Use monolithic velocity range for color scaling (partitioned may differ)
vmax_mono = max(np.max(v[1]) for v in mono_vel)

# Create velocity animation
print("  Rendering velocity video...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

def plot_frame(frame_idx):
    for ax in axes:
        ax.clear()

    pts_m, vmag_m = mono_vel[frame_idx]
    pts_p, vmag_p = part_vel[frame_idx]

    # Each subplot uses its own normalization for best visibility
    sc1 = axes[0].scatter(pts_m[:, 2], pts_m[:, 0], c=vmag_m, s=2,
                          vmin=0, vmax=vmax_mono, cmap='coolwarm')
    axes[0].set_title(f'Monolithic FSI')
    axes[0].set_xlabel('z')
    axes[0].set_ylabel('x')
    axes[0].set_aspect('equal')

    sc2 = axes[1].scatter(pts_p[:, 2], pts_p[:, 0], c=vmag_p, s=2,
                          vmin=0, vmax=vmax_mono, cmap='coolwarm')
    axes[1].set_title(f'Partitioned FSI')
    axes[1].set_xlabel('z')
    axes[1].set_ylabel('x')
    axes[1].set_aspect('equal')

    fig.suptitle(f'Velocity Magnitude (step {frame_idx+1}/{n_steps}, dt=1e-4)', fontsize=14)
    return sc1, sc2

plot_frame(0)
plt.tight_layout()

ani = animation.FuncAnimation(fig, plot_frame, frames=n_steps, interval=200, blit=False)
ani.save(os.path.join(out_dir, 'velocity_comparison.gif'), writer='pillow', fps=5)
plt.close()
print(f"  Saved {out_dir}/velocity_comparison.gif")

# Create pressure comparison at final step
print("  Creating final pressure comparison...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

mono_final = meshio.read(os.path.join(mono_dir, f"result_{n_steps:03d}.vtu"))
part_final = meshio.read(os.path.join(part_dir, f"result_{n_steps:03d}.vtu"))

p_mono = mono_final.point_data.get('Pressure', np.zeros(len(mono_final.points)))
p_part = part_final.point_data.get('Pressure', np.zeros(len(part_final.points)))
pmax = max(np.max(np.abs(p_mono)), np.max(np.abs(p_part)))

sc1 = axes[0].scatter(mono_final.points[:, 2], mono_final.points[:, 0],
                       c=p_mono, s=2, vmin=-pmax, vmax=pmax, cmap='RdBu_r')
axes[0].set_title('Monolithic FSI')
axes[0].set_xlabel('z')
axes[0].set_ylabel('x')
axes[0].set_aspect('equal')
plt.colorbar(sc1, ax=axes[0], label='Pressure')

sc2 = axes[1].scatter(part_final.points[:, 2], part_final.points[:, 0],
                       c=p_part, s=2, vmin=-pmax, vmax=pmax, cmap='RdBu_r')
axes[1].set_title('Partitioned FSI')
axes[1].set_xlabel('z')
axes[1].set_ylabel('x')
axes[1].set_aspect('equal')
plt.colorbar(sc2, ax=axes[1], label='Pressure')

fig.suptitle(f'Pressure at step {n_steps}', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'pressure_comparison.png'), dpi=150)
plt.close()
print(f"  Saved {out_dir}/pressure_comparison.png")

print("\nDone!")
