#!/bin/bash

#SBATCH --job-name=collect-sim-data
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#
#######################################


bash velocity_vs_activity.sh
