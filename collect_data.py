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

if __name__ == "__main__":
    collect_zeta_data()

