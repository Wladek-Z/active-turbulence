#!/bin/bash

# Paths
PARAM_FILE="in_turbulence275000"
OUTPUT_DIR="../velocity vs time 64"
SIM_BINARY="./Qrho275000"      # Simulation executable
OUTPUT_FILE="velocity275000.dat"     # Simulation output file

# Create results directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Loop over X parameter sets
for i in $(seq 1 20); do

    # Step 1: Edit parameter file    
    if [[ $i -le 40 ]]; then
    	activity=$(echo "$i * 0.0001"| bc -l)
    else
    	activity=$(echo "$i * 0.0002 - 0.004" | bc -l)
    fi
    	
    parameter_value="$activity, 0"
    echo "$parameter_value" > "$PARAM_FILE"

    # Step 2: Start simulation
    $SIM_BINARY < "$PARAM_FILE"

    # Step 3: Rename and move output
    NEW_NAME="velocity_${activity}.dat"
    if [[ -f "$OUTPUT_FILE" ]]; then
        mv "$OUTPUT_FILE" "$OUTPUT_DIR/$NEW_NAME"
    else
        echo "Warning: Output file not found for run $i"
    fi
done

