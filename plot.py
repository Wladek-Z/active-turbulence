import numpy as np
import matplotlib.pyplot as plt

system = "64"

from pathlib import Path
local_path = Path(__file__).parent / f"velocity vs time {system}"

def plot_vz(data):
    """Plot average velocity vs activity parameter"""
    z = data[:, 0]
    v = data[:, 1]

    xmin = np.min(z)
    xmax = np.max(z)
    ymin = np.min(v)
    ymax = np.max(v)

    plt.figure(figsize=(8, 6))
    plt.scatter(z, v, marker='+',  color='b', s=10)
    #plt.plot(z, v, linestyle='-', color='b')
    plt.xlabel(r'activity parameter ($\zeta$)')
    plt.ylabel(r'average velocity [$su$]')
    plt.title(rf'Velocity vs Activity (${system} \times {system}$ system)')
    #plt.xlim(0.002, 0.0045) 
    #plt.ylim(0, 0.0025)
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


if __name__ == "__main__":
    file_path = local_path / 'SD_data.txt'
    data = np.loadtxt(file_path)
    data = data[::5]
    plot_SD(data)