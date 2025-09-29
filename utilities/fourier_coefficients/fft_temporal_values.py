import numpy as np
from matplotlib import pyplot as plt
import os 
pi = np.pi

def fft(nt, temporal_values, d, n):
    """
    Computes the Fourier Coefficients for a given set of temporal values. Line by line translation of the fft function in svMultiphysics (fft.cpp)

    Parameters:
        nt (int): Total number of time points in temporal_values file
        temporal_values (array): List of time and corresponding data values size. Each row of array is [time, d1, d2, d3]
        d (int): Dimensions of the data (usually 1 for temporal values)
        n (int): Number of Fourier components to be computed

    Returns:
        dict: A dictionary containing the Fourier coefficients and other related values:
            - 'r': Real part of the Fourier coefficients
            - 'i': Imaginary part of the Fourier coefficients
            - 'qi': Initial values for each data dimension
            - 'qs': Slopes for each data dimension
            - 'ti': Initial time
            - 'T': Total duration of the time series
    """
    t = np.zeros(nt)
    q = np.zeros((d, nt))

    # Extract time and data values
    for i in range(nt):                         # for each time point (row)
        t[i] = temporal_values[i][0]            # save the time value
        for j in range(d):                      # for each dimension (column)  
            q[j, i] = temporal_values[i][j + 1] # save the corresponding data value

    ti = t[0]              # initial time
    T = t[-1] - t[0]       # total duration

    # Extract the initial values and the linear slopes of the data
    qi = q[:, 0].copy()
    qs = (q[:, -1] - q[:, 0]) / T

    # Pre-processing: de-trending the data by removing the initial value and linear slope
    for i in range(nt):
        t[i] -= ti
        for j in range(d):
            q[j, i] -= qi[j] + qs[j] * t[i]

    # Initialize real and imaginary output arrays of the fourier coefficients 
    r = np.zeros((d, n))
    i_ = np.zeros((d, n))

    for n_idx in range(n):
        tmp = float(n_idx)

        for i in range(nt - 1):
            ko = 2.0 * pi * tmp * t[i] / T
            kn = 2.0 * pi * tmp * t[i + 1] / T

            for j in range(d):
                s = (q[j, i + 1] - q[j, i]) / (t[i + 1] - t[i])     # local slope
                if n_idx == 0:          # For DC component, use trapezoidal rule
                    r[j, n_idx] += 0.5 * (t[i + 1] - t[i]) * (q[j, i + 1] + q[j, i])
                else:                   
                    r[j, n_idx] += s * (np.cos(kn) - np.cos(ko))
                    i_[j, n_idx] -= s * (np.sin(kn) - np.sin(ko))

        if n_idx == 0:                  # For DC component, scale by T  
            r[:, n_idx] /= T
        else:                           # For other components, scale by T/(pi^2 * tmp^2)
            scale = 0.5 * T / (pi * pi * tmp * tmp)
            r[:, n_idx] *= scale
            i_[:, n_idx] *= scale

    return {
        'r': r,
        'i': i_,
        'qi': qi,
        'qs': qs,
        'ti': ti,
        'T': T
    }

def write_fourier_coeff_file(filename, result, d, n):
    """
    Writes the Fourier coefficients and related information to the .fcs file format.

    Parameters:
        filename (str): The name of the file to write to.
        result (dict): The result dictionary containing Fourier coefficients and related values.
        d (int): The number of spatial dimensions.
        n (int): The number of Fourier components.
    """

    with open(filename, 'w') as f:
        # Write initial time and total duration
        f.write(f"{result['ti']} {result['T']}\n")

        # Write qi and qs for each dimension
        for j in range(d):
            f.write(f"{result['qi'][j]} {result['qs'][j]}\n")

        # Write number of Fourier components
        f.write(f"{n}\n")

        # Write real and imaginary parts
        for k in range(n):
            real_parts = ' '.join(f"{result['r'][j][k]}" for j in range(d))
            imag_parts = ' '.join(f"{result['i'][j][k]}" for j in range(d))
            # format the output to be 16 decimal places in scientific notation
            # also add 4 spaces between real and imaginary parts
            real_parts = ' '.join(f"{float(part):.16e}" for part in real_parts.split())
            imag_parts = ' '.join(f"{float(part):.16e}" for part in imag_parts.split())
            f.write(f"{real_parts}    {imag_parts}\n")
    return -1 

def recon_fft(qi, qs, ti, T, r, i, times, nfcs):
    """
    Reconstructs the signal using Fourier coefficients. Works currently only for one dimension

    Parameters:
        qi (float): Initial values for each data dimension
        qs (float): Slopes for each data dimension
        ti (float): Initial time
        T (float): Total duration of the time series
        r (array): Real part of the Fourier coefficients
        i (array): Imaginary part of the Fourier coefficients
        times (array): Array of time points at which to reconstruct the signal
        nfcs (int): Number of Fourier components

    Returns:
        array: Reconstructed signal values at the specified time points
    """
    t_rel = np.array(times - ti)
    rec = float(qi)*np.ones(len(t_rel)) + float(qs) * t_rel + r[0]*np.ones(len(t_rel))  # DC + linear trend
    for k in range(1, nfcs):                                    # Skip k=0 since already added
        freq = 2 * np.pi * k * t_rel / T
        fn = r[k] * np.cos(freq) - i[k] * np.sin(freq)
        rec += fn
    return rec 

def visualize_fft(result, temporal_values, nfcs):
    """
    Visualizes the original and reconstructed signals using FFT.

    Parameters:
        result (dict): The result dictionary containing Fourier coefficients and related values.
        temporal_values (array): The original temporal values.
        nfcs (int): The number of Fourier components.
    """
    times = temporal_values[:, 0]
    ti, T, qi, qs, r, i = result['ti'], result['T'], result['qi'][0], result['qs'][0], result['r'][0], result['i'][0]
    rec = recon_fft(qi, qs, ti, T, r, i, times, nfcs)

    plt.scatter(times, temporal_values[:, 1], label='Temporal Values')
    plt.plot(times, rec, label='Interpolated fcs', linestyle='--', color='red')
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Amplitude', fontsize=12)
    plt.title('FFT Reconstruction', fontsize=14)
    plt.legend(fontsize=12)
    plt.savefig("fft_reconstruction.png", dpi=400)
    plt.show()

if __name__ == "__main__":
    # Suppose temporal_values is a list of [time, val1, val2, ...]
    file_values = np.loadtxt(os.path.join('./lumen_inlet.flow'))

    # Extract the first row to get number of timepoints and number of fourier components
    nt, nf = file_values[0]
    nt, nf = int(nt), int(nf)
    # Extract the time and temporal values 
    temporal_values = file_values[1:nt+1, :]
    d = len(temporal_values[0]) - 1             # dimensions: size of each row minus 1 for the time data

    result = fft(nt, temporal_values, d, nf)

    print("Real part (r):", result['r'])
    print("Imag part (i):", result['i'])
    print("Initial values (qi):", result['qi'])
    print("Slopes (qs):", result['qs'])
    print("Initial time (ti):", result['ti'])
    print("Total duration (T):", result['T'])

    # After computing `result` with your fft function:
    output_filename = os.path.join('./lumen_inlet.fcs')
    write_fourier_coeff_file(output_filename, result, d, nf)

    print(f"Fourier coefficient file written to {output_filename}")
    visualize_fft(result, temporal_values, nf)
