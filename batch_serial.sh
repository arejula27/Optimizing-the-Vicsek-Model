#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Second job allocation
#SBATCH -J serial
#SBATCH -t 1:00:00
#SBATCH -A edu24.DD2356
#SBATCH -p main 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -e error_file_serial.e
#SBATCH --nodes=1

output_file=serial_output.txt

> $output_file
for i in $(seq 1 5); do
    echo "Run $i - Mode: serial" >> $output_file
    srun bin/bench.out serial >> $output_file
    echo "" >> $output_file    # add blank line
done