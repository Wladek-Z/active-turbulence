#!/bin/bash

#SBATCH --job-name=collect-sim-data-32
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#
#######################################


bash velocity_vs_activity.sh
