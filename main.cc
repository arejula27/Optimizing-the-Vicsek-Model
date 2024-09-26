/**
 * @file main.cc
 * @brief Run the flocking bird simulation and plot the results
 */

#include <iostream>
#include <random>
#include <vector>

#include "simulation.h"

#ifdef USE_MPI_FUNCTION
// Define the flag USE_MPI_FUNCTION to use the MPI optimized function
#define simulation_function simulation_mpi
#elif defined(USE_OPENMP_FUNCTION)
// Define the flag USE_OPENMP_FUNCTION to use the OpenMP optimized function
#define simulation_function simulation_openmp
#else
// If USE_MPI_FUNCTION and USE_OPENMP_FUNCTION are not defined, the standard function will be used
#define simulation_function simulation
#endif



/**
 * @brief Main function for running the flocking bird simulations and plotting the results.
 * 
 * @return Integer indicating the exit status of the program.
 */
int main()
{
    // Finite Volume simulation

    // Simulation parameters
    double v0 = 1.0;  // Initial velocity
    double eta = 0.5; // random fluctuation in angle (in radians)
    double L = 10;    // size of the box
    double R = 1;     // interaction radius
    double dt = 0.2;  // time step
    int Nt = 500;     // number of time steps
    int N = 500;      // number of birds

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

    // Simulation
    for (int i = 0; i < Nt; i++)
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

        theta = simulation_function(x, y, theta, R);

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

    // plot
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe)
    {
        // Configure the picture
        fprintf(gnuplotPipe, "set terminal png\n");
        fprintf(gnuplotPipe, "set output 'assets/birds.png'\n");
        fprintf(gnuplotPipe, "set xrange [0:%f]\n", L);
        fprintf(gnuplotPipe, "set yrange [0:%f]\n", L);
        fprintf(gnuplotPipe, "unset border\n");
        fprintf(gnuplotPipe, "unset xtics\n");
        fprintf(gnuplotPipe, "unset ytics\n");
        fprintf(gnuplotPipe, "set size square\n");

        // Draw  vectors
        for (int i = 0; i < N; ++i)
        {
            fprintf(gnuplotPipe, "set arrow %d from %f,%f to %f,%f lt 1 lw 0 filled\n", i + 1, x[i], y[i], x[i] + vx[i], y[i] + vy[i]);
        }
        fprintf(gnuplotPipe, "plot NaN title ''\n"); // Esto es necesario para que se muestren las flechas
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
        std::cout << "Pictutre saved as'birds.png'" << std::endl;
    }
    else
    {
        std::cerr << "Error: Gnuplot pipe could not be opened." << std::endl;
        return 1;
    }
    return 0;
}
