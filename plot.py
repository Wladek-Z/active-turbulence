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
    plt.plot(z, SD, linestyle='-', color='r', label=f'{system}x{system} system')
    plt.fill_between(z, SD - SD_err, SD + SD_err, color='r', alpha=0.2)
    plt.xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    plt.ylabel(r'$\sigma (\nu)$ [$su$]', fontsize=12)
    #plt.title(rf'Standard Deviation of Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.show()

def plot_vt(data):
    """Plot average velocity against timestep"""
    t = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(t, v, linestyle='-', color='g', label=f'{system}x{system} system')
    plt.xlabel('timestep', fontsize=12)
    plt.ylabel(r'$\overline{v}$ [$su$]', fontsize=12)
    #plt.title(rf'Average Velocity Evolution (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.show()

def plot_mean_v(data):
    """Plot mean velocity in steady-state against activity parameter"""
    z = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(z, v, linestyle='-', color='b', label=f'{system}x{system} system')
    plt.xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    plt.ylabel(r'$\langle v \rangle$ in steady-state [$su$]', fontsize=12)
    #plt.title(rf'Mean Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend(loc='lower right')
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

def plot_both_mean_vz(data32, data64):
    """Plot mean velocity in steady-state vs activity parameter for both 32x32 and 64x64 systems."""
    data32 = np.array(data32)
    data64 = np.array(data64)

    z32 = data32[:, 0]
    v32 = data32[:, 1]

    z64 = data64[:, 0]
    v64 = data64[:, 1]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Mean Velocity in Steady-State vs. Activity', fontsize=16)

    ax[0].plot(z32, v32, label='32x32 system', color='b')
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'$\langle v \rangle$ in steady-state [su]', fontsize=12)
    ax[0].legend(loc='lower right')
    ax[0].grid(True)

    ax[1].plot(z64, v64, label='64x64 system', color='b')
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'mean velocity in steady-state [$su$]', fontsize=12)
    ax[1].legend(loc='lower right')
    ax[1].grid(True)

    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

def plot_both_SD(data32, data64):
    """Plot standard deviation of velocity vs activity parameter for both 32x32 and 64x64 systems."""
    data32 = np.array(data32)
    data64 = np.array(data64)

    z32 = data32[:, 0]
    SD32 = data32[:, 1]
    SD_err32 = data32[:, 2]

    z64 = data64[:, 0]
    SD64 = data64[:, 1]
    SD_err64 = data64[:, 2]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Standard Deviation of Velocity vs. Activity', fontsize=16)

    ax[0].plot(z32, SD32, label='32x32 system', color='r')
    ax[0].fill_between(z32, SD32 - SD_err32, SD32 + SD_err32, color='r', alpha=0.2)
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'$\sigma(v)$ [$su$]', fontsize=12)
    ax[0].legend(loc='lower right')
    ax[0].grid(True)

    ax[1].plot(z64, SD64, label='64x64 system', color='r')
    ax[1].fill_between(z64, SD64 - SD_err64, SD64 + SD_err64, color='r', alpha=0.2)
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'standard deviation of velocity [$su$]', fontsize=12)
    ax[1].legend(loc='lower right')
    ax[1].grid(True)
    
    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

if __name__ == "__main__":
    """
    system = "32"
    local_path = Path(__file__).parent / f"velocity vs time AB0.0 {system}"
    file_path = local_path / 'velocity_.0900.dat'

    data = np.loadtxt(file_path)[::1]
    
    plot_vt(data)
    """
    local_path32 = Path(__file__).parent / "velocity vs time AB0.0 32"
    local_path64 = Path(__file__).parent / "velocity vs time AB0.0 64"

    file_path32 = local_path32 / 'mean_v_data.txt'
    file_path64 = local_path64 / 'mean_v_data.txt'

    data32 = np.loadtxt(file_path32)[::4]
    data64 = np.loadtxt(file_path64)[::2]

    plot_both_mean_vz(data32, data64)
    