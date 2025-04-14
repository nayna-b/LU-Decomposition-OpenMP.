#!/bin/bash
#SBATCH --export=ALL
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --time=00:29:00 
#SBATCH --partition=interactive
#SBATCH --reservation=comp422
#SBATCH --job-name=lu-scaling
#SBATCH --output=lu-scaling.out
#SBATCH --error=lu-scaling.err

export KMP_AFFINITY=verbose # just to print the bindings

MATRIX_SIZE=8000
OUTPUT_FILE=result/scaling_times.txt

mkdir -p result

for ((i = 1; i <= 32; i++)); do 
  echo "----- Threads = $i -----" >> $OUTPUT_FILE
  OMP_NUM_THREADS=$i OMP_PROC_BIND=spread ./lu-omp $MATRIX_SIZE $i >> $OUTPUT_FILE
  echo "" >> $OUTPUT_FILE
  echo "Completed run with $i threads"
done
