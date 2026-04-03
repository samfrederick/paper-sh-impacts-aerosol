#!/bin/bash

# Note - enable execution first ('chmod +x compile_nsh.sh') prior to running this script

# source your conda.sh pointed to by the conda base environment
# https://github.com/conda/conda/issues/7980#issuecomment-441358406
#source /Users/sf20/opt/anaconda3/etc/profile.d/conda.sh
source ~/anaconda3/etc/profile.d/conda.sh

conda init bash
conda activate paper_env
f2py -c nsh.f90 -m nsh