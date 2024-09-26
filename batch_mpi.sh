#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Second job allocation
#SBATCH -J MPI_4
#SBATCH -t 4:00:00
#SBATCH -A edu24.DD2356
#SBATCH -p main 
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -e error_file_MPI.e
#SBATCH --nodes=4

nodes=4     # SET THIS EQUAL TO NODES FLAG VALUE

# Function to run the executable multiple times
run_executable() {
    local map="mpi_data"
    local mode=$1
    local num_iterations=$2
    local output_file=$3
    local proc_count=$4

    # create directory if it does not exist
    if [ ! -d "${map}" ]; then
        mkdir "${map}"
    fi

    # remove old values
    > "${map}/${output_file}"

    for i in $(seq 1 $num_iterations); do
        echo "Run $i - Mode: $mode" >> "${map}/${output_file}"
        srun -n $proc_count bin/bench.out $mode >> "${map}/${output_file}"
        echo "" >> "${map}/${output_file}"    # add blank line
    done
}

for processes in 1 2 4 8 16 32 64 128; do
    echo "Running mpi with $processes processes:"
    # Run with mode "MPI"
    run_executable "MPI" 5 "${nodes}_nodes_${processes}_processes.txt" $processes
done



