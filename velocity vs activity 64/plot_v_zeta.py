import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
local_path = Path(__file__).parent / "velocity vs activity 64"

def load_data(file_path):
    """Load data from a text file."""
    data = np.loadtxt(file_path)
    return data

def plot_v_z(data):
    """Plot average velocity vs activity parameter"""
    zeta = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.scatter(zeta, v, marker='x',  color='b', s=10)
    plt.plot(zeta, v, linestyle='--', color='b', alpha=0.5)
    plt.xlabel('Activity parameter (zeta)')
    plt.ylabel(r'Average velocity [$ms^{-1}$]')
    plt.title('Velocity vs activity')
    #plt.xlim(0, 0.008) # for 64x64 system
    #plt.ylim(-0.001, 0.015) # for 64x64 system
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = local_path / 'v_vs_zeta_data.txt'
    data = load_data(file_path)
    plot_v_z(data)