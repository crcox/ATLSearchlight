#!/bin/bash

# -- SGE optoins (whose lines must begin with #$)

#$ -S /bin/bash # The jobscript is written for the bash shell
#$ -V # Inherit environment settings (e.g., from loaded modulefiles)
#$ -cwd # Run the job in the current directory
#$ -l short

# Relevant modules? Perhaps?
# compilers/gcc/4.6.2
# compilers/gcc/4.7.0 * most current that is officially compatible with mex
# compilers/gcc/4.8.2
# compilers/gcc/4.9.0
# compilers/gcc/6.3.0
#
# apps/binapps/matlab/R2010a
# apps/binapps/matlab/R2011a
# apps/binapps/matlab/R2012a
# apps/binapps/matlab/R2013a
# apps/binapps/matlab/R2014a
# apps/binapps/matlab/R2015a
# apps/binapps/matlab/R2015aSP1
# -- the commands to be executed (programs to be run) on a compute node:

module load compilers/gcc/4.7.0
module load apps/binapps/matlab/R2015aSP1
# `which matlab` will return path/bin/matlab, so 2x dirname will strip off
# bin/matlab and leave the path.
export MATLABDIR="$(dirname $(dirname $(which matlab)))"
make
