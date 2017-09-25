#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M jlstiles@berkeley.edu 
#$ -m beas
export OMP_NUM_THREADS=1

module load gcc/4.9.3

# Load Java if any R packages need RJava (bartMachine, h2o, etc.)
module load java

# Load a better linear algebra system.
module load lapack

# GPU computation modules. CUDA is 7.5, cudnn is 4.0.
module load cuda cudnn


