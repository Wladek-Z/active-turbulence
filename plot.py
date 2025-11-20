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
    plt.ylabel(r'average velocity [su]')
    plt.title(rf'Velocity vs Activity (${system} \times {system}$)')
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
    plt.ylabel(r'$v_{rms}$ in steady-state [su]')
    plt.title(rf'Root Mean Square Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.show()

def plot_SD(data):
    """Plot standard deviation of velocity vs activity parameter with shaded error region"""
    data1 = data[::4]  # downsample for clarity
    z = data1[:, 0]
    SD = data1[:, 1]
    SD_err = data1[:, 2]

    plt.figure(figsize=(8, 6))
    plt.plot(z, SD, linestyle='-', color='r', label=rf'{system}$\times${system}')
    plt.fill_between(z, SD - SD_err, SD + SD_err, color='r', alpha=0.2)
    plt.xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    plt.ylabel(r'$\sigma (v)$ [su]', fontsize=12)
    #plt.title(rf'Standard Deviation of Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.show()

def plot_vt(data):
    """Plot average velocity against timestep"""
    t = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(t, v, linestyle='-', color='g', label=rf'${system} \times {system}$')
    plt.xlabel('timestep', fontsize=12)
    plt.ylabel(r'$\overline{v}$ [su]', fontsize=12)
    #plt.title(rf'Average Velocity Evolution (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_mean_v(data):
    """Plot mean velocity in steady-state against activity parameter"""
    data1 = data[::2]
    z = data1[:, 0]
    v = data1[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(z, v, linestyle='-', color='b', label=rf'${system} \times {system}$')
    plt.xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    plt.ylabel(r'$\langle v \rangle$ in steady-state [su]', fontsize=12)
    #plt.title(rf'Mean Velocity vs Activity (${system} \times {system}$ system)')
    plt.grid(True)
    plt.legend()
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
    data32 = np.array(data32)[::4]
    data64 = np.array(data64)[::2]

    z32 = data32[:, 0]
    v32 = data32[:, 1]

    z64 = data64[:, 0]
    v64 = data64[:, 1]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Mean Velocity in Steady-State vs. Activity', fontsize=16)

    ax[0].plot(z32, v32, label=r'32$\times$32', color='b')
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'$\langle v \rangle$ in steady-state [su]', fontsize=12)
    ax[0].set_ylabel(r'$\langle v \rangle$ in steady-state [su]', fontsize=12)
    ax[0].legend()
    ax[0].grid(True)

    ax[1].plot(z64, v64, label=r'64$\times$64', color='b')
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'mean velocity in steady-state [$su$]', fontsize=12)
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

def plot_both_SD(data32, data64):
    """Plot standard deviation of velocity vs activity parameter for both 32x32 and 64x64 systems."""
    data32 = np.array(data32)[::8]
    data64 = np.array(data64)[::4]

    z32 = data32[:, 0]
    SD32 = data32[:, 1]
    SD_err32 = data32[:, 2]

    z64 = data64[:, 0]
    SD64 = data64[:, 1]
    SD_err64 = data64[:, 2]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Standard Deviation of Velocity vs. Activity', fontsize=16)

    ax[0].plot(z32, SD32, label=r'32$\times$32', color='r')
    ax[0].fill_between(z32, SD32 - SD_err32, SD32 + SD_err32, color='r', alpha=0.2)
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'$\sigma(v)$ [su]', fontsize=12)
    ax[0].legend(loc='upper left')
    ax[0].grid(True)

    ax[1].plot(z64, SD64, label=r'64$\times$64', color='r')
    ax[1].fill_between(z64, SD64 - SD_err64, SD64 + SD_err64, color='r', alpha=0.2)
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'standard deviation of velocity [$su$]', fontsize=12)
    ax[1].legend(loc='upper left')
    ax[1].grid(True)
    
    plt.tight_layout()  # Add extra padding at bottom
    plt.show()


def plot_both_error(data32, data64):
    """Plot standard deviation of velocity vs activity parameter for both 32x32 and 64x64 systems."""
    data32 = np.array(data32)[::8]
    data64 = np.array(data64)[::4]

    z32 = data32[:, 0]
    SD32 = data32[:, 1]
    SD_err32 = data32[:, 2]

    z64 = data64[:, 0]
    SD64 = data64[:, 1]
    SD_err64 = data64[:, 2]

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    #fig.suptitle(rf'Standard Deviation of Velocity vs. Activity', fontsize=16)

    ax[0].plot(z32, SD_err32, color='r')
    ax[0].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    ax[0].set_ylabel(r'error on $\sigma(v)$ [su]', fontsize=12)
    ax[0].legend(loc='upper left')
    ax[0].grid(True)

    ax[1].plot(z64, SD_err64, color='r')
    ax[1].set_xlabel(r'activity parameter ($\zeta$)', fontsize=12)
    #ax[1].set_ylabel(r'standard deviation of velocity [$su$]', fontsize=12)
    ax[1].legend(loc='upper left')
    ax[1].grid(True)
    
    plt.tight_layout()  # Add extra padding at bottom
    plt.show()

def plot_both_vt_fft(data1, data2):
    """plot the velocity vs time data in the top row, and the FFT in the bottom row for two datasets."""
    t1 = data1[:, 0]
    v1 = data1[:, 1]

    t2 = data2[:, 0]
    v2 = data2[:, 1]

    # Compute the original FFT for both datasets
    dft1 = rfft(v1[-50:])
    dft1[0] = 0
    dft1 = np.abs(dft1)
    freq1 = rfftfreq(len(t1[-50:]), d=(t1[1] - t1[0]))

    dft2 = rfft(v2[-50:])
    dft2[0] = 0
    dft2 = np.abs(dft2)
    freq2 = rfftfreq(len(t2[-50:]), d=(t2[1] - t2[0]))

    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    # Top row: velocity vs time
    ax[0, 0].plot(t1, v1, linestyle='-', color='g')
    ax[0, 0].set_xlabel('timestep', fontsize=12)
    ax[0, 0].set_ylabel(r'$\overline{v}$ [su]', fontsize=12)
    ax[0, 0].set_title(r"$\zeta = 0.0200$")
    ax[0, 0].grid(True)

    ax[0, 1].plot(t2, v2, linestyle='-', color='g')
    ax[0, 1].set_xlabel('timestep', fontsize=12)
    ax[0, 1].set_title(r"$\zeta = 0.0800$")
    ax[0, 1].grid(True)

    # Bottom row: FFT plots
    ax[1, 0].plot(freq1, dft1, linestyle='-', color='orange')
    ax[1, 0].set_xlabel(r'Frequency [cycles timestep$^{-1}$]')
    ax[1, 0].set_ylabel('|DFT|')
    ax[1, 0].grid(True)

    ax[1, 1].plot(freq2, dft2, linestyle='-', color='orange')
    ax[1, 1].set_xlabel(r'Frequency [cycles timestep$^{-1}$]')
    ax[1, 1].grid(True)

    plt.tight_layout()  # Add extra padding at bottom
    plt.show()


if __name__ == "__main__":
    
    system = "64"
    local_path = Path(__file__).parent / f"velocity vs time AB0.1 {system}"
    file_path = local_path / 'velocity_.0040.dat'

    data = np.loadtxt(file_path)
    
    plot_vt(data)
    """
    local_path32 = Path(__file__).parent / "velocity vs time AB0.0 32"
    local_path64 = Path(__file__).parent / "velocity vs time AB0.0 64"

    file_path32 = local_path32 / 'SD_data.txt'
    file_path64 = local_path64 / 'SD_data.txt'

    data32 = np.loadtxt(file_path32)
    data64 = np.loadtxt(file_path64)

    plot_both_error(data32, data64)
    """"""
    system = "32"
    local_path = Path(__file__).parent / f"velocity vs time AB0.0 {system}"
    file_path1 = local_path / 'velocity_.0200.dat'
    file_path2 = local_path / 'velocity_.0800.dat'

    data1 = np.loadtxt(file_path1)
    data2 = np.loadtxt(file_path2)

    plot_both_vt_fft(data1, data2)
    """