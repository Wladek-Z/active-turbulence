#!/bin/bash

#SBATCH --job-name=vt-64-275000
#SBATCH --partition=long
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#
#######################################


bash v_vs_t275000.sh
