#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M jlstiles@berkeley.edu 
#$ -m beas
#$ OMP_NUM_THREADS=1

R --vanilla < longdelta3_10000.R > longdelta3_10000.Rout
