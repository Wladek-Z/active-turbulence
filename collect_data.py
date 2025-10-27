from pathlib import Path
#local_path = Path(__file__).parent / "velocity vs time 32"

import numpy as np
import re
from scipy.stats import bootstrap

def collect_zeta_data():
    data = []  # will contain lists of [zeta, average speed]

    for velocity in local_path.rglob("velocity_*.dat"):
        # Extract the float number from the filename
        filename = velocity.name  # e.g., "velocity_.002.dat"
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)

        if match:
            zeta = float(match.group(1))

            with open(velocity, 'r') as file:
                lines = file.readlines()
                parts = lines[-1].strip().split()
                speed = float(parts[1])
            data.append([zeta, speed])
        else:
            print(f"Warning: No valid data found in file {velocity.name}")

    # Save the collected data to a text file
    with open(local_path / "vz_data.txt", 'w') as outfile:
        outfile.write("# zeta    average_speed\n")
        for zeta, avg_speed in sorted(data):
            outfile.write(f"{zeta} {avg_speed}\n")

def collect_SD_data():
    """Read last 200 lines of data, calculate the standard deviation of velocity for each file. Save to SD_data.txt"""
    data = []  # will contain lists of [zeta, SD, SD_err]

    for velocity in local_path.rglob("velocity_*.dat"):
        # Extract the float number from the filename
        filename = velocity.name  # e.g., "velocity_.002.dat"
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)
        K = 200  # Number of bootstrap resamples

        if match:
            zeta = float(match.group(1))

            with open(velocity, 'r') as file:
                lines = file.readlines()[-200:]  # Read last 200 lines
                speeds = [float(line.strip().split()[1]) for line in lines]
                SD = np.std(speeds)
                res = bootstrap(
                    (speeds,),
                    statistic=std_func,
                    n_resamples=K,
                    method='percentile',
                    random_state=42
                )
                SD_err = res.standard_error

            data.append([zeta, SD, SD_err])

        else:
            print(f"Warning: No valid data found in file {velocity.name}")

    # Save the collected data to a text file
    with open(local_path / "SD_data.txt", 'w') as outfile:
        outfile.write("# zeta    SD    SD_err\n")
        for zeta, SD, SD_err in sorted(data):
            outfile.write(f"{zeta} {SD} {SD_err}\n")

def std_func(x):
    return np.std(x, ddof=1)

def find_Kmax(data, threshold=0.01, max_resamples=100000):
    """Find number of bootstrap resamples (Kmax) to achieve relative error < threshold."""
    speeds = data[:, 1][-200:]
    orig_std = np.std(speeds, ddof=1)

    K = 100  # start with 100 resamples
    step = 100
    last_se = None

    while K <= max_resamples:
        res = bootstrap(
            (speeds,),
            statistic=std_func,
            n_resamples=K,
            method='percentile',
            random_state=42
        )
        
        se = res.standard_error

        if last_se is not None:
            rel_error = abs(se - last_se) / se
            if rel_error < threshold:
                print(f"Converged with K={K}, std error={se}, original std={orig_std}")
                return K
        
        last_se = se
        K += step

    print(f"Did not converge after {max_resamples} resamples")

def mean_v_data():
    """Read last 200 lines of data, calculate the mean velocity for each file. Save to mean_v_data.txt"""
    data = []  # will contain lists of [zeta, mean_speed]

    for velocity in local_path.rglob("velocity_*.dat"):
        # Extract the float number from the filename
        filename = velocity.name  # e.g., "velocity_.002.dat"
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)

        if match:
            zeta = float(match.group(1))

            with open(velocity, 'r') as file:
                lines = file.readlines()[-200:]  # Read last 200 lines
                speeds = [float(line.strip().split()[1]) for line in lines]
                mean_speed = np.mean(speeds)

            data.append([zeta, mean_speed])
            
        else:
            print(f"Warning: No valid data found in file {velocity.name}")

    # Save the collected data to a text file
    with open(local_path / "mean_v_data.txt", 'w') as outfile:
        outfile.write("# zeta    mean_speed\n")
        for zeta, mean_speed in sorted(data):
            outfile.write(f"{zeta} {mean_speed}\n")

if __name__ == "__main__":
    local_path = Path(__file__).parent / "velocity vs time AB0.1 64"
    collect_zeta_data()
    mean_v_data()
    collect_SD_data()
    