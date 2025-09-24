import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
local_path = Path(__file__).parent

def load_data(file_path):
    """Load data from a text file."""
    data = np.loadtxt(file_path)
    return data

def plot_v_z(data):
    """Plot average velocity vs activity parameter"""
    zeta = data[:, 0]
    v = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(zeta, v, marker='o', linestyle='-', color='b')
    plt.xlabel('Activity parameter (zeta)')
    plt.ylabel(r'Average velocity [$ms^{-1}$]')
    plt.title('Velocity vs activity')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = local_path / 'v_vs_zeta_data.txt'
    data = load_data(file_path)
    plot_v_z(data)