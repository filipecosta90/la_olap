#!/bin/bash

# Preparing environment for Parallel Studio XE 2016 and enabling ifort 16
echo "Loading Parallel Studio XE 2016 Compilers ..."
module purge
cd 
cd /share/apps/intel/parallel_studio_xe_2016.0.047/compilers_and_libraries_2016/linux/bin
source compilervars_global.sh intel64
cd /share/apps/intel/parallel_studio_xe_2016.0.047/compilers_and_libraries_2016/linux/mkl/bin
source mklvars.sh intel64
echo "DONE!"

echo "Ready now for running make!"
