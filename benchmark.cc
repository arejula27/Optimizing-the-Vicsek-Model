/**
 * @file benchmark.cc
 * @brief Benchmarks the elapsed time for running the flocking bird simulation with the different theta implementations
 */

#include <iostream>
#include <assert.h>
#include <random>
#include <chrono>
#include <functional>

#include "simulation.h" 

const int Nt = 500; // number of time steps
const int N = 2000;  // number of birds

/**
 * @brief Perform benchmark simulations using a specified simulation function.
 * 
 * @param rank The rank of the process, 0 if only 1 process is used.
 * @param sim_function The simulation function to be used for the benchmark.
 * @param func_string A string describing the simulation function used.
 * @return A vector containing the final theta values after the simulation.
 */
std::vector<double> benchmark_simulation(int rank, std::function<std::vector<double>(const std::vector<double>&, const std::vector<double>&, std::vector<double>, double)> sim_function, std::string func_string)
{
    // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step

    // Initialise
    std::default_random_engine generator;
    generator.seed(static_cast<unsigned int>(17));
    std::uniform_real_distribution<double> uniform(0, 1);

    // Initialize bird position
    std::vector<double> x(N);
    std::vector<double> y(N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = uniform(generator) * L;
        y[i] = uniform(generator) * L;
    }

    // Initialize bird velocity
    std::vector<double> theta(N);
    std::vector<double> vx(N);
    std::vector<double> vy(N);
    for (int i = 0; i < N; ++i)
    {
        theta[i] = 2 * M_PI * uniform(generator);
        vx[i] = v0 * cos(theta[i]);
        vy[i] = v0 * sin(theta[i]);
    }

    // start time counter
    auto start = std::chrono::high_resolution_clock::now();

    // Simulation
    for (int j = 0; j < Nt; j++)
    {
        for (int i = 0; i < N; ++i)
        {
            // Update position
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;

            // Apply periodic boundary conditions
            x[i] = fmod(x[i], L);
            y[i] = fmod(y[i], L);
            if (x[i] < 0)
                x[i] += L; // Handle negative values
            if (y[i] < 0)
                y[i] += L; // Handle negative values
        }

        theta = sim_function(x, y, theta, R);

        // Update velocities
        for (int i = 0; i < N; ++i)
        {
            // Update theta randomly
            theta[i] += eta * (uniform(generator) - 0.5);
            // Update velocity components
            vx[i] = v0 * cos(theta[i]);
            vy[i] = v0 * sin(theta[i]);
        }
    }
    // stop time counter

    if(rank == 0){
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << func_string << std::endl;
        std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    }

    return theta;
}

/**
 * @brief Main function for running benchmark simulations.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return Integer indicating the exit status of the program.
 */
int main(int argc, char* argv[])
{
    int provided, rank;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get each process rank

    if(argc > 1)
    {
        if(std::string(argv[1]) == "serial")
        {
            benchmark_simulation(rank, simulation, "simulation");
        } 
        else if(std::string(argv[1]) == "openMP")
        {
            benchmark_simulation(rank, simulation_openmp, "simulation_openmp");
        } 
        else if(std::string(argv[1]) == "openMP_dy")
        {
            benchmark_simulation(rank, simulation_openmp_dy, "simulation_openmp_dy");
        }
        else if(std::string(argv[1]) == "MPI")
        {
            benchmark_simulation(rank, simulation_mpi, "simulation_mpi");
        } else {
            std::cout << "Parameter must be either serial, openMP, openMP_dy or MPI" << std::endl;
        }
    } else { //Run all benchmarks
        std::vector<double> theta1(N);
        std::vector<double> theta2(N);
        if(rank == 0){
            std::cout << "Running all benchmarks..." << std::endl;

            theta1 = benchmark_simulation(rank, simulation, "simulation");

            theta2 = benchmark_simulation(rank, simulation_openmp, "simulation_openmp");
            assert(theta1.size() == theta2.size());
            assert(theta1 == theta2);

            theta2 = benchmark_simulation(rank, simulation_openmp_dy, "simulation_openmp_dy");
            assert(theta1.size() == theta2.size());
            assert(theta1 == theta2);
        }

        MPI_Barrier(MPI_COMM_WORLD); //Wait for all processes before running MPI simulations

        theta2 = benchmark_simulation(rank, simulation_mpi, "simulation_mpi");

        if(rank == 0){
            assert(theta1.size() == theta2.size());
            assert(theta1 == theta2);       
        }

        if(rank == 0)
            std::cout << "All simulations gave the same theta" << std::endl;
    }

    MPI_Finalize();
}