#!/bin/bash

#SBATCH --job-name=collect-sim-data-32
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#
#######################################


bash velocity_vs_activity.sh
