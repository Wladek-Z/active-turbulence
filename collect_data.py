from pathlib import Path
local_path = Path(__file__).parent / "velocity vs time 64"

import numpy as np
import re

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

def collect_std_data():
    """Read last 200 lines of data, calculate the standard deviation of velocity for each file. Save to std_data.txt"""
    data = []  # will contain lists of [zeta, std_dev]

    for velocity in local_path.rglob("velocity_*.dat"):
        # Extract the float number from the filename
        filename = velocity.name  # e.g., "velocity_.002.dat"
        match = re.search(r"velocity_([-+]?\d*\.?\d+)\.dat", filename)

        if match:
            zeta = float(match.group(1))

            with open(velocity, 'r') as file:
                lines = file.readlines()[-200:]  # Read last 200 lines
                speeds = [float(line.strip().split()[1]) for line in lines]
                std_dev = np.std(speeds)
            data.append([zeta, std_dev])
        else:
            print(f"Warning: No valid data found in file {velocity.name}")

    # Save the collected data to a text file
    with open(local_path / "std_data.txt", 'w') as outfile:
        outfile.write("# zeta    std_dev\n")
        for zeta, std_dev in sorted(data):
            outfile.write(f"{zeta} {std_dev}\n")

if __name__ == "__main__":
    local_path = Path(__file__).parent / "velocity vs time 32"
    #collect_zeta_data()
    collect_std_data()
