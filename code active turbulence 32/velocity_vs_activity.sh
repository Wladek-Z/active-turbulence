#!/bin/bash

# Paths (adjust as needed)
PARAM_FILE="in_turbulence"
OUTPUT_DIR="../velocity vs activity 32"
SIM_BINARY="./Qrho"      # Your simulation executable
OUTPUT_FILE="Qtensor.txt"     # File your simulation produces

# Create results directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Loop over X parameter sets
for i in $(seq 1 250); do

    # Step 1: Edit parameter file
    activity=$(echo "$i * 0.0002"| bc -l)
    parameter_value="$activity, 0"
    echo "$parameter_value" > "$PARAM_FILE"

    # Step 2: Start simulation with variable timeout
    if [[ $i -le 49 ]]; then
        timeout 300s $SIM_BINARY < "$PARAM_FILE"
    else
        timeout 60s $SIM_BINARY < "$PARAM_FILE"
    fi
    SIM_PID=$!

    # Step 3: Rename and move output
    NEW_NAME="Qtensor_${activity}.txt"
    if [[ -f "$OUTPUT_FILE" ]]; then
        mv "$OUTPUT_FILE" "$OUTPUT_DIR/$NEW_NAME"
    else
        echo "Warning: Output file not found for run $i"
    fi
done
