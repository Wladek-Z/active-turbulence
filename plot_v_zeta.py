import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
local_path = Path(__file__).parent / "velocity vs activity 32"

def load_data(file_path):
    """Load data from a text file."""
    data = np.loadtxt(file_path)
    return data

def plot_v_z(data):
    """Plot average velocity vs activity parameter"""
    zeta = data[:, 0]
    v = data[:, 1]

    xmin = np.min(zeta)
    xmax = np.max(zeta)
    ymin = np.min(v)
    ymax = np.max(v)

    plt.figure(figsize=(8, 6))
    #plt.scatter(zeta, v, marker='+',  color='b', s=10)
    plt.plot(zeta, v, linestyle='-', color='b')
    plt.xlabel('Activity parameter (zeta)')
    plt.ylabel(r'Average velocity [$ms^{-1}$]')
    plt.title(r'Velocity vs activity ($32 \times 32$)')
    plt.xlim(xmin, 0.007) 
    plt.ylim(ymin, 0.006)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = local_path / 'v_vs_zeta_data.txt'
    data = load_data(file_path)
    plot_v_z(data)