import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
local_path = Path(__file__).parent / "velocity vs time 32"

def load_data(file_path):
    """Load data from a text file."""
    data = np.loadtxt(file_path)
    return data

def plot_v_z(data):
    """Plot average velocity vs activity parameter"""
    t = data[:, 0]
    v = data[:, 1]

    xmin = np.min(t)
    xmax = np.max(t)
    ymin = np.min(v)
    ymax = np.max(v)

    plt.figure(figsize=(8, 6))
    #plt.scatter(zeta, v, marker='+',  color='b', s=10)
    plt.plot(t, v, linestyle='-', color='b')
    plt.xlabel(r'Time [$su$]')
    plt.ylabel(r'Average velocity [$su$]')
    plt.title(r'Velocity vs Time ($32 \times 32$)')
    #plt.xlim(xmin, 0.004) 
    #plt.ylim(0, 0.01)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = local_path / 'velocity_.0500.dat'
    data = load_data(file_path)
    plot_v_z(data)