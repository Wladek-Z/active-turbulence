#!/bin/bash

#SBATCH --job-name=collect-vt-data-32
#SBATCH --partition=long
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#
#######################################


bash velocity_vs_time.sh
