import numpy as np
import matplotlib.pyplot as plt

system = "32"

from pathlib import Path
local_path = Path(__file__).parent / f"velocity vs time {system}"

def load_data(file_path):
    """Load data from a text file."""
    data = np.loadtxt(file_path)
    return data

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
    plt.xlim(0.002, 0.0045) 
    plt.ylim(0, 0.0025)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = local_path / 'vz_data.txt'
    data = load_data(file_path)
    plot_vz(data)