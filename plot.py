import numpy as np
from scipy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
from pathlib import Path


def plot_vz(data):
    """Plot average velocity vs activity parameter"""
    z = data[:, 0]
    v = data[:, 1]

    xmin = np.min(z)
    xmax = np.max(z)
    ymin = np.min(v)
    ymax = np.max(v)

    plt.figure(figsize=(8, 6))
    #plt.scatter(z, v, marker='+',  color='b', s=10)
    plt.plot(z, v, linestyle='-', color='b')
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'average velocity [$su$]')
    plt.title(rf'Velocity vs Activity (${system} \times {system}$ system)')
    plt.xlim(0.0002, 0.007) 
    plt.ylim(0, 0.006)
    plt.grid(True)
    plt.show()

def plot_v_rms(data):
    """Plot root mean square velocity vs activity parameter"""
    z = data[:, 0]
    v = data[:, 1]

    xmin = np.min(z)
    xmax = np.max(z)
    ymin = np.min(v)
    ymax = np.max(v)

    plt.figure(figsize=(8, 6))
    plt.plot(z, v, linestyle='-', color='b')
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'$v_{rms}$ in steady-state [$su$]')
    plt.title(rf'Root Mean Square Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

def plot_SD(data):
    """Plot standard deviation of velocity vs activity parameter with shaded error region"""
    z = data[:, 0]
    SD = data[:, 1]
    SD_err = data[:, 2]

    plt.figure(figsize=(8, 6))
    plt.plot(z, SD, linestyle='-', color='r')
    plt.fill_between(z, SD - SD_err, SD + SD_err, color='r', alpha=0.2)
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'standard deviation of velocity [$su$]')
    plt.title(rf'Standard Deviation of Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

def plot_vt(data):
    """Plot average velocity against timestep"""
    t = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(t, v, linestyle='-', color='g')
    plt.xlabel('timestep')
    plt.ylabel(r'average velocity [$su$]')
    plt.title(rf'Average Velocity Evolution (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

def plot_mean_v(data):
    """Plot mean velocity in steady-state against activity parameter"""
    z = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(z, v, linestyle='-', color='b')
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'mean velocity in steady-state [$su$]')
    plt.title(rf'Mean Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

def plot_fft(data):
    """Plot FFT of velocity data"""
    d = data[-50:]
    t = d[:, 0]
    v = d[:, 1]

    # Compute the original FFT
    dft = rfft(v)
    
    # Manually set first component to zero
    dft[0] = 0
    dft = np.abs(dft)

    freq = rfftfreq(len(t), d=(t[1] - t[0]))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Discrete Fourier Transform of Velocity (${system} \times {system}$ system)', fontsize=16)

    ax[0].plot(freq, dft, linestyle='-', color='orange')
    ax[0].set_xlabel(r'Frequency [cycles timestep$^{-1}$]')
    ax[0].set_ylabel('|DFT|')
    #ax[0].set_title(rf'Discrete Fourier Transform of Velocity (${system} \times {system}$ system)')
    ax[0].grid(True)

    ax[1].semilogy(freq, dft, linestyle='-', color='orange')
    ax[1].set_xlabel(r'Frequency [cycles timestep$^{-1}$]')
    ax[1].set_ylabel('|DFT| (log scale)')
    #ax[1].set_title(rf'Discrete Fourier Transform of Velocity (${system} \times {system}$ system)')
    ax[1].grid(True)
    
    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

if __name__ == "__main__":
    system = "32"
    local_path = Path(__file__).parent / f"velocity vs time AB0.0 {system}"
    file_path = local_path / 'SD_data.txt'

    data = np.loadtxt(file_path)[::8]
    
    plot_SD(data)