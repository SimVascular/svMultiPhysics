"""
This script plots the displacement and pressure-volume loop for the LV_HolzapfelOgden_active case 
and compares results to those found in the benchmark paper by Aróstica et al. (2025).

The resulting plots are included in this directory. If the user wishes to generate the plots,
they can do so by downloading the benchmark dataset from https://zenodo.org/records/14260459
and running solver.xml with 1000 time steps to obtain the svMultiPhysics results.
"""

import pickle
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pyvista as pv
from scipy.interpolate import interp1d

sns.set_theme()

DATA_PATH = Path(__file__).parent / "data"

# if DATA_PATH doesn't exist, print error messsage and exit
if not DATA_PATH.exists():
    print(f"Error: Benchmark path {DATA_PATH} does not exist. Benchmark results can be found at https://zenodo.org/records/14260459")
    exit(1)

# Path to svMultiPhysics results
results_path = 'results'
if not Path(results_path).exists():
    print(f"Error: svMultiPhysics results path {results_path} does not exist. Run solver.xml with 1000 time steps to generate results.")
    exit(1)


def load_dataset_from_pickle(filename: str | Path) -> dict[str, np.ndarray]:
    if not isinstance(filename, Path):
        filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError(f"File {filename} not found")

    with open(filename.as_posix(), "rb") as fl:
        return pickle.load(fl)


def get_array_at_point(mesh, point, array_name):
    """Interpolate array value at an arbitrary point using PyVista's sample method."""
    # Create a PolyData object for the query point
    point_poly = pv.PolyData(np.array([point]))
    # Interpolate the array at the point
    sampled = point_poly.sample(mesh)
    if array_name in sampled.point_data:
        return sampled.point_data[array_name][0]
    else:
        raise ValueError(f"Array '{array_name}' not found in mesh")


def get_displacements_at_points(start_timestep, end_timestep, step, timestep_size, results_folder, sample_points):
    """Get displacements at specified points over time."""
    t = []
    displacements = []
    
    for k in range(start_timestep, end_timestep + 1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))
        
        # Get displacement at each sample point
        point_displacements = []
        for point in sample_points:
            if k == 0:
                # At time 0, displacement is zero
                disp = np.array([0.0, 0.0, 0.0])
            else:
                # Get displacement at the point
                disp = get_array_at_point(result, point, 'Displacement')
            point_displacements.append(disp)
        
        t.append(k * timestep_size)
        displacements.append(point_displacements)
    
    return np.array(t), np.array(displacements)


def calc_volume_3D(start_timestep, end_timestep, step, timestep_size, results_folder, reference_surface, save_intermediate_data=False, intermediate_output_folder=None):
    """Calculate volume over time using 3D mesh."""
    t = []
    vol = []
    
    # Load reference surface
    ref_lumen = pv.read(reference_surface)
    
    for k in range(start_timestep, end_timestep + 1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))
        
        # Sample result onto ref_lumen
        resampled_lumen = ref_lumen.sample(result)
        
        # Warp resampled surface by displacement
        warped_lumen = resampled_lumen.warp_by_vector('Displacement')
        
        # Save warped and filled lumen if requested
        if save_intermediate_data:
            warped_lumen.save(os.path.join(intermediate_output_folder, f'resampled_warped_and_filled_{k:03d}.vtp'))
        
        # Add time and volume to arrays
        t.append(k * timestep_size)
        vol.append(warped_lumen.volume)
        
        print(f"Iteration: {k}, Volume: {warped_lumen.volume}")
    
    return (t, vol)


def get_svmultiphysics_displacement_data():
    """Get svMultiPhysics displacement data in the same format as benchmark datasets."""
    # Set points needed to calculate the displacements
    p0 = np.array([0.025, 0.03, 0.0])
    p1 = np.array([0.0, 0.03, 0.0])
    sample_points = [p0, p1]
    results_path = 'results'
    
    # Get the displacements at the sample points
    t_displacements, displacements = get_displacements_at_points(10, 1000, 10, 1e-3, results_path, sample_points)
    
    # Format data to match benchmark dataset structure
    svmultiphysics_data = {
        "time": t_displacements,
        "displacement": {
            "p0": {
                "ux": displacements[:, 0, 0],
                "uy": displacements[:, 0, 1], 
                "uz": displacements[:, 0, 2]
            },
            "p1": {
                "ux": displacements[:, 1, 0],
                "uy": displacements[:, 1, 1],
                "uz": displacements[:, 1, 2]
            }
        }
    }
    
    return svmultiphysics_data


def plot_pv_loop_only():
    """Plot only the pressure-volume loop."""
    # Plot P-V loop
    results_path = 'results'
    t, vol = calc_volume_3D(10, 1000, 10, 1e-3, results_path, 'mesh/mesh-surfaces/endo.vtp')
    vol = np.array(vol) * 1e6  # Convert to mL
    
    # Load pressure data from the second column of pressure.dat
    pressure_data = np.loadtxt('pressure.dat')
    pressure_time = pressure_data[1:, 0]  # First column: time
    pressure_values = pressure_data[1:, 1] / 133.322387415  # Second column: pressure
    
    # Interpolate pressure to match volume time points
    pressure_interp = interp1d(pressure_time, pressure_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    pressure_at_volume_times = pressure_interp(t)
    
    plt.figure(figsize=(10, 8))
    plt.plot(vol, pressure_at_volume_times, 'b-', linewidth=2, label='Pressure-Volume Loop')
    plt.xlabel('Volume [mL]', fontsize=12)
    plt.ylabel('Pressure [mmHg]', fontsize=12)
    plt.title('Pressure-Volume Loop', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('p-v_loop.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print some statistics
    print(f"Maximum volume: {np.max(vol):.2f} mL")
    print(f"Minimum volume: {np.min(vol):.2f} mL")
    print(f"Stroke volume: {np.max(vol) - np.min(vol):.2f} mL")
    print(f"Maximum pressure: {np.max(pressure_at_volume_times):.2f} mmHg")
    print(f"Minimum pressure: {np.min(pressure_at_volume_times):.2f} mmHg")


def plot_displacements_and_pv_loop():
    """Plot displacements and pressure-volume loop from process_results.py."""
    # Get svMultiPhysics data
    svmultiphysics_data = get_svmultiphysics_displacement_data()
    
    # Load benchmark data
    data = np.load('Step_1_US_P1_h5.npz')
    coords = ['x', 'y', 'z']
    
    # Plot the displacements for both points and for each of three coordinates on 2 x 3 grid
    plt.figure(figsize=(12, 12))
    for i in range(3):
        plt.subplot(3, 2, 2*i+1)
        plt.plot(svmultiphysics_data['time'], svmultiphysics_data['displacement']['p0'][f'u{coords[i]}'], label='svMultiPhysics')
        plt.plot(data['time'], data['u_0'][:, i], label='benchmark')
        plt.xlabel('Time [s]')
        plt.ylabel(coords[i] + ' Displacement [m]')
        plt.legend()
        plt.tight_layout()
        
        plt.subplot(3, 2, 2*i+2)
        plt.plot(svmultiphysics_data['time'], svmultiphysics_data['displacement']['p1'][f'u{coords[i]}'], label='svMultiPhysics')
        plt.plot(data['time'], data['u_1'][:, i], label='benchmark')
        plt.xlabel('Time [s]')
        plt.ylabel(coords[i] + ' Displacement [m]')
        plt.legend()
        plt.tight_layout()
    
    # Label column 1 'p_0' and column 2 'p_1'
    plt.subplot(3, 2, 1)
    plt.title(r'$p_0$')
    plt.subplot(3, 2, 2)
    plt.title(r'$p_1$')
    
    plt.tight_layout()
    plt.savefig('displacements.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Plot P-V loop
    results_path = 'results'
    t, vol = calc_volume_3D(0, 1000, 10, 1e-3, results_path, 'mesh/mesh-surfaces/endo.vtp')
    vol = np.array(vol) * 1e6  # Convert to mL
    
    # Load pressure data from the second column of pressure.dat
    pressure_data = np.loadtxt('pressure.dat')
    pressure_time = pressure_data[1:, 0]  # First column: time
    pressure_values = pressure_data[1:, 1] / 133.322387415  # Second column: pressure
    
    # Interpolate pressure to match volume time points
    pressure_interp = interp1d(pressure_time, pressure_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    pressure_at_volume_times = pressure_interp(t)
    
    plt.figure(figsize=(10, 8))
    plt.plot(vol, pressure_at_volume_times, 'b-', linewidth=2, label='Pressure-Volume Loop')
    plt.xlabel('Volume [mL]', fontsize=12)
    plt.ylabel('Pressure [mmHg]', fontsize=12)
    plt.title('Pressure-Volume Loop', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('p-v_loop.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print some statistics
    print(f"Maximum volume: {np.max(vol):.2f} mL")
    print(f"Minimum volume: {np.min(vol):.2f} mL")
    print(f"Stroke volume: {np.max(vol) - np.min(vol):.2f} mL")
    print(f"Maximum pressure: {np.max(pressure_at_volume_times):.2f} mmHg")
    print(f"Minimum pressure: {np.min(pressure_at_volume_times):.2f} mmHg")


LABEL_NAMES = [
    "CARPentry",
    "Ambit",
    "4C",
    "Simula",
    "CHimeRA",
    r"$\mathcal{C}$Heart",
    r"life$^{\mathbf{X}}$",
    r"SimVascular $\mathbb{P}_1$",
    r"SimVascular $\mathbb{P}_2$",
    "COMSOL",
]

LABEL_NAMES_BIV = [
    "CARPentry",
    "Ambit",
    "4C",
    "Simula",
    "CHimeRA",
    r"$\mathcal{C}$Heart",
    r"life$^{\mathbf{X}}$",
    r"SimVascular $\mathbb{P}_1$",
    "COMSOL",
]

COLORS = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:cyan",
    "tab:olive",
]

COLORS_BIV = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
]

LABELS_P1 = [
    r"Simula-$\mathbb{P}_1$",
    r"CHimeRA-$\mathbb{P}_1$",
    r"$\mathcal{C}$Heart-$\mathbb{P}_1$",
    r"life$^{\mathbf{X}}$-$\mathbb{P}_1$",
    r"COMSOL-$\mathbb{P}_1$",
]

LABELS_P2 = [
    r"Simula-$\mathbb{P}_2$",
    r"CHimeRA-$\mathbb{P}_2$",
    r"$\mathcal{C}$Heart-$\mathbb{P}_2$",
    r"life$^{\mathbf{X}}$-$\mathbb{P}_2$",
    r"COMSOL-$\mathbb{P}_2$",
]

COLORS_P1_P2 = ["tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive"]


TEAMS_DATASETS_1 = [
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_carpentry.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_ambit.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_4c.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_simula.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_chimera.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_cheart.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_lifex.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_simvascular_p1p1.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_simvascular_p2.pickle"
    ),
    load_dataset_from_pickle(
        DATA_PATH / "monoventricular_nonblinded_step_1_group_comsol.pickle"
    ),
]


def compute_plot_displacement_monoventricular(
    teams_datasets: list[dict[str, np.ndarray]],
    filename: str | Path,
    labels_names: list[str],
    colors: list[str],
) -> None:
    """Computes the displacement figures based on the groups datasets."""
    if isinstance(filename, str):
        filename = Path(filename)

    if not isinstance(teams_datasets, list):
        raise Exception("teams_datasets must be a list")

    if not len(teams_datasets) == len(labels_names):
        raise Exception("teams_datasets and labels_names must have the same length")

    if not len(teams_datasets) == len(colors):
        raise Exception("teams_datasets and colors must have the same length")

    filename.parent.mkdir(parents=True, exist_ok=True)

    fig, axs = plt.subplots(3, 2, sharex=True, figsize=(12, 12))

    axs[0, 0].set_title(r"Particle $\mathbf{p}_0$")
    axs[0, 0].set_ylabel("Displacement x-component [m]")
    axs[1, 0].set_ylabel("Displacement y-component [m]")
    axs[2, 0].set_ylabel("Displacement z-component [m]")
    axs[2, 0].set_xlabel("Time [s]")

    axs[0, 1].set_title(r"Particle $\mathbf{p}_1$")
    axs[2, 1].set_xlabel("Time [s]")

    for data, lbl, color in zip(teams_datasets, labels_names, colors):
        for i, u_type in zip([0, 1, 2], ["ux", "uy", "uz"]):
            if lbl == "svMultiPhysics":
                axs[i, 0].plot(
                    data["time"], data["displacement"]["p0"][u_type], label=lbl, color="k", linestyle="--"
                )
                axs[i, 1].plot(
                    data["time"], data["displacement"]["p1"][u_type], label=lbl, color="k", linestyle="--"
                )
            else:
                axs[i, 0].plot(
                    data["time"], data["displacement"]["p0"][u_type], label=lbl, color=color
                )
                axs[i, 1].plot(
                    data["time"], data["displacement"]["p1"][u_type], label=lbl, color=color
                )

    plt.legend(loc="best", fancybox=True, shadow=True)
    fig.tight_layout()
    fig.savefig(filename.as_posix(), bbox_inches="tight", dpi=120)
    # plt.show()


def compute_statistics(
    datasets: list[dict[str, np.ndarray]], point: str = "p0"
) -> dict[str, np.ndarray]:
    time = datasets[0]["time"]
    number_of_element_per_dataset = datasets[0]["time"].shape[0]
    number_of_datasets = len(datasets)
    condensated_dataset_ux = np.zeros(
        (number_of_element_per_dataset, number_of_datasets)
    )
    condensated_dataset_uy = np.zeros_like(condensated_dataset_ux)
    condensated_dataset_uz = np.zeros_like(condensated_dataset_ux)
    red_dataset_u = np.zeros((number_of_datasets,))

    for i, dataset in enumerate(datasets):
        condensated_dataset_ux[:, i] = dataset["displacement"][point]["ux"]
        condensated_dataset_uy[:, i] = dataset["displacement"][point]["uy"]
        condensated_dataset_uz[:, i] = dataset["displacement"][point]["uz"]

    mean_ux = np.mean(condensated_dataset_ux, axis=1)
    mean_uy = np.mean(condensated_dataset_uy, axis=1)
    mean_uz = np.mean(condensated_dataset_uz, axis=1)

    for i in range(number_of_datasets):
        diff_to_mean_norm = np.sqrt(
            np.abs(condensated_dataset_ux[:, i] - mean_ux) ** 2
            + np.abs(condensated_dataset_uy[:, i] - mean_uy) ** 2
            + np.abs(condensated_dataset_uz[:, i] - mean_uz) ** 2
        )

        mean_norm = np.sqrt(
            np.abs(mean_ux) ** 2 + np.abs(mean_uy) ** 2 + np.abs(mean_uz) ** 2
        )

        red_dataset_u[i] = np.mean(diff_to_mean_norm / mean_norm, axis=0)

    statistics_dataset = {
        "time": time,
        "mean_ux": mean_ux,
        "mean_uy": mean_uy,
        "mean_uz": mean_uz,
        "std_ux": np.std(condensated_dataset_ux, axis=1),
        "std_uy": np.std(condensated_dataset_uy, axis=1),
        "std_uz": np.std(condensated_dataset_uz, axis=1),
        "red_u": red_dataset_u,
    }

    return statistics_dataset


def compute_plots_monoventricular_nonblinded_step_1() -> None:
    """Computes displacement curves for the monoventricular non-blinded step 1"""
    # Get svMultiPhysics data
    svmultiphysics_data = get_svmultiphysics_displacement_data()
    
    # Add svMultiPhysics to the datasets, labels, and colors
    all_datasets = TEAMS_DATASETS_1 + [svmultiphysics_data]
    all_labels = LABEL_NAMES + ["svMultiPhysics"]
    all_colors = COLORS + ["tab:red"]  # Add red color for svMultiPhysics
    
    compute_plot_displacement_monoventricular(
        all_datasets,
        "./comparison_plots_p0_p1_step_1_nonblinded.png",
        all_labels,
        all_colors,
    )


def compute_statistics_per_dataset(
    labels: list[str], datasets: list[dict[str, np.ndarray]], header: str = ""
) -> None:
    if len(labels) != len(datasets):
        raise ValueError("labels and datasets do not match in length")

    print(header)
    statistics_p0 = compute_statistics(datasets, point="p0")
    statistics_p1 = compute_statistics(datasets, point="p1")

    red_p0 = np.round(statistics_p0["red_u"], decimals=3)
    red_p1 = np.round(statistics_p1["red_u"], decimals=3)

    print(f"Team {'TEAM NAME':25s} - {'p0':4s}, {'p1':4s}")
    for label, i in zip(labels, range(len(red_p0))):
        print(f"Team {label:25s} - {red_p0[i]:4.3f}, {red_p1[i]:4.3f}")

    print("-----------------------------------\n")


if __name__ == "__main__":
    compute_plots_monoventricular_nonblinded_step_1()
    
    # Also plot PV loop separately
    print("Plotting pressure-volume loop...")
    plot_pv_loop_only()
