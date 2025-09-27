from pathlib import Path
local_path = Path(__file__).parent / "velocity vs activity 32"

import numpy as np
import re

data = [] # will contain lists of [zeta, average speed]

for qtensor in local_path.rglob("Qtensor_*.txt"):
    # Extract the float number from the filename
    filename = qtensor.name  # e.g., "Qtensor_.002.txt"
    match = re.search(r"Qtensor_([-+]?\d*\.?\d+)\.txt", filename)

    if match:
        zeta = float(match.group(1))
        speeds = []

        with open(qtensor, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()  # Remove whitespace/newlines
                if line:  # Skip empty lines
                    values = list(map(float, line.split()))
                    speeds.append(np.linalg.norm(values[8:11]))  # Assuming velocity components are at indices 8, 9, 10

        data.append([zeta, np.mean(speeds)])
        
with open(local_path / "v_vs_zeta_data.txt", 'w') as outfile:
    outfile.write("# zeta    average_speed\n")
    for zeta, avg_speed in sorted(data):
        outfile.write(f"{zeta} {avg_speed}\n")


