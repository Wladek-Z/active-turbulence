#!/bin/bash

#SBATCH --job-name=vt-64-175000
#SBATCH --partition=long
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#
#######################################


bash v_vs_t175000.sh
