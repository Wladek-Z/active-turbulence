import scipy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re

def get_data(directory):
    """Collect velocity and timestep data for each zeta value from files in the specified directory.
       If spacing > 1, take every 'spacing' data point to reduce data size."""
    zeta = []
    velocity_list = []
    timestep_list = []

    for vel_file in directory.rglob("velocity_*.dat"):
        # Extract the float number from the filename
        filename = vel_file.name  # e.g., "velocity_.002.dat"
        #print(f"Processing file: {filename}")
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)

        if match:
            zeta_val = float(match.group(1))

            with open(vel_file, 'r') as file:
                lines = file.readlines()[-50:]  # Read last 50 lines
                timesteps = [float(line.strip().split()[0]) for line in lines]
                velocities = [float(line.strip().split()[1]) for line in lines]

            zeta.append(zeta_val)
            velocity_list.append(velocities)
            timestep_list.append(timesteps)

        else:
            print(f"Warning: No valid data found in file {vel_file.name}")
    
    order = np.argsort(zeta)
    zeta = np.array(zeta)[order]
    velocity_list = [velocity_list[i] for i in order]
    timestep_list = [timestep_list[i] for i in order]

    return zeta, timestep_list, velocity_list

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
        FT[0] = 0  # Remove zero-frequency component
        FT = np.abs(FT)
        nu = fft.rfftfreq(N, dt)

        # Find average frequency <nu>, average square frequency <nu^2>
        sum_FT_squared = np.sum(FT**2)
        
        if sum_FT_squared == 0:
            print(f"Warning: Sum of FT^2 is {sum_FT_squared} for zeta={zeta}. Signal likely too constant.")
            avg_nu = avg_nu2 = 0
            #start = i+1 # start from next index
        else:
            avg_nu = np.sum(nu * FT**2) / sum_FT_squared
            avg_nu2 = np.sum(nu**2 * FT**2) / sum_FT_squared

        # Find standard deviation of frequency spectrum
        SD_nu = np.sqrt(avg_nu2 - avg_nu**2)

        SD_freq.append([zeta, SD_nu])
    
    return SD_freq#[start:] # [0]: zeta, [1]: SD of frequency

def plot_SD_freq(SD_freq, spacing=1):
    """Plot standard deviation of frequency vs activity parameter"""
    SD_freq = np.array(SD_freq)
    z = SD_freq[:, 0]
    SD_nu = SD_freq[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(z[::spacing], SD_nu[::spacing], color='m', label=f'{system}x{system} system')
    plt.xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    plt.ylabel(r'$\sigma(\nu)$ [cycles timestep$^{-1}$]', fontsize=12)
    #plt.title(rf'Standard Deviation of Frequency vs. Activity', fontsize=16)
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.show()

def plot_both(data32, data64, spacing32, spacing64):
    """Plot standard deviation of frequency vs activity parameter for both 32x32 and 64x64 systems."""
    data32 = np.array(data32)
    data64 = np.array(data64)

    z32 = data32[:, 0]
    SD_nu32 = data32[:, 1]

    z64 = data64[:, 0]
    SD_nu64 = data64[:, 1]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Standard Deviation of Frequency vs. Activity', fontsize=16)

    ax[0].plot(z32[::spacing32], SD_nu32[::spacing32], label='32x32 system', color='m')
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'$\sigma(\nu)$ [cycles timestep$^{-1}$]', fontsize=12)
    ax[0].legend(loc='lower right')
    ax[0].grid(True)

    ax[1].plot(z64[::spacing64], SD_nu64[::spacing64], label='64x64 system', color='m')
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'standard deviation of frequency [cycles timestep$^{-1}$]')
    ax[1].legend(loc='lower right')
    ax[1].grid(True)
    
    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

if __name__ == "__main__":
    
    system = "64"
    local_path = Path(__file__).parent / f"velocity vs time AB0.1 {system}"
    
    if system == "32":
        spacing = 4
    elif system == "64":
        spacing = 2

    z, t_list, v_list = get_data(local_path)
    SD_freq = frequency_analysis(z, t_list, v_list)
    plot_SD_freq(SD_freq, spacing=spacing)
    """
    local_path32 = Path(__file__).parent / "velocity vs time AB0.0 32"
    local_path64 = Path(__file__).parent / "velocity vs time AB0.0 64"

    z32, t_list32, v_list32 = get_data(local_path32)
    SD_freq32 = frequency_analysis(z32, t_list32, v_list32)

    z64, t_list64, v_list64 = get_data(local_path64)
    SD_freq64 = frequency_analysis(z64, t_list64, v_list64)

    plot_both(SD_freq32, SD_freq64, spacing32=4, spacing64=2)
    """
