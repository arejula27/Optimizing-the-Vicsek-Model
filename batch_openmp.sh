#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Second job allocation
#SBATCH -J openMP
#SBATCH -t 4:00:00
#SBATCH -A edu24.DD2356
#SBATCH -p main 
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --nodes=1
#SBATCH -e error_file_openMP.e

# Function to run the executable multiple times
run_executable() {
    local mode=$1
    local num_iterations=$2
    local output_file=$3
    local map="${mode}_data"

    # create directory if it does not exist
    if [ ! -d "${map}" ]; then
        mkdir "${map}"
    fi

    # remove old values
    > "${map}/${output_file}"

    for i in $(seq 1 $num_iterations); do
        echo "Run $i - Mode: $mode" >> "${map}/${output_file}"
        srun bin/bench.out $mode >> "${map}/${output_file}"
        echo "" >> "${map}/${output_file}"    # add blank line
    done
}

# Outer loop over number of threads = powers of 2
for mode in "openMP" "openMP_dy"; do
    for num_threads in 1 2 4 8 16 32 64 128; do
        echo "Running $mode with $num_threads threads:"
        export OMP_NUM_THREADS=$num_threads
        run_executable $mode 5 "${num_threads}.txt"
    done
done