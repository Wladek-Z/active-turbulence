import scipy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re

def get_data(directory, spacing=1):
    """Collect velocity and timestep data for each zeta value from files in the specified directory.
       If spacing > 1, take every 'spacing' data point to reduce data size."""
    zeta = []
    velocity_list = []
    timestep_list = []
    i = -1

    for vel_file in directory.rglob("velocity_*.dat"):
        # Check if the current index is a multiple of the spacing
        # Note: this works equivalently to data[::spacing]
        i += 1
        if i % spacing != 0:
            continue

        # Extract the float number from the filename
        filename = vel_file.name  # e.g., "velocity_.002.dat"
        #print(f"Processing file: {filename}")
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)

        if match:
            zeta_val = float(match.group(1))

            with open(vel_file, 'r') as file:
                lines = file.readlines()[-200:]  # Read last 200 lines
                timesteps = [float(line.strip().split()[0]) for line in lines]
                velocities = [float(line.strip().split()[1]) for line in lines]

            zeta.append(zeta_val)
            velocity_list.append(velocities)
            timestep_list.append(timesteps)

        else:
            print(f"Warning: No valid data found in file {vel_file.name}")

    return np.array(zeta), np.array(timestep_list), np.array(velocity_list)

def frequency_analysis(z, t_list, v_list):
    """Perform Fourier analysis to find the standard deviation of the frequency spectrum for each zeta value."""
    SD_freq = []
    
    for i, zeta in enumerate(z):
        t = t_list[i]
        v = v_list[i]

        dt = t[1] - t[0]
        N = len(v)

        # Perform Fourier Transform
        FT = fft.rfft(v)
        nu = fft.rfftfreq(N, dt)

        # Find average frequency <nu>, average square frequency <nu^2>
        avg_nu = np.sum(nu * np.abs(FT)**2) / np.sum(np.abs(FT)**2)
        avg_nu2 = np.sum(nu**2 * np.abs(FT)**2) / np.sum(np.abs(FT)**2)

        # Find standard deviation of frequency spectrum
        SD_nu = np.sqrt(avg_nu2 - avg_nu**2)

        SD_freq.append([zeta, SD_nu])

    return SD_freq # [0]: zeta, [1]: SD of frequency

def plot_SD_freq(SD_freq):
    """Plot standard deviation of frequency vs activity parameter"""
    SD_freq = np.array(SD_freq)
    z = SD_freq[:, 0]
    SD_nu = SD_freq[:, 1]
    
    order = np.argsort(z)
    z = z[order]
    SD_nu = SD_nu[order]

    plt.figure(figsize=(8, 6))
    plt.plot(z, SD_nu, color='m')
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'standard deviation of frequency [cycles timestep$^{-1}$]')
    plt.title(rf'Standard Deviation of Frequency vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    system = "64"
    local_path = Path(__file__).parent / f"velocity vs time {system}"
    
    if system == "32":
        spacing = 10
    elif system == "64":
        spacing = 5

    z, t_list, v_list = get_data(local_path, spacing=spacing)
    SD_freq = frequency_analysis(z, t_list, v_list)
    plot_SD_freq(SD_freq)
