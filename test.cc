/**
 * @file test.cc
 * @brief Executes some tests to verify the correctness of the simulation functions
 */

#include <iostream>
#include <assert.h>

#include "simulation.h"

const double R = 1.0;
const double threshold = 0.0001;

void init_vectors(std::vector<double> &x, std::vector<double> &y, std::vector<double> &theta)
{
    x = {1.0, 0.0, 1.0, 1.0};
    y = {1.0, 1.0, 0.0, 0.3};
    theta = {1.0, 1.0, 1.0, 1.0};
}

struct test_case
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> theta;
    std::vector<double> correct_result;
};

test_case test_cases[] = {
    {{1.0, 0.0, 1.0, 1.0}, {1.0, 1.0, 0.0, 0.3}, {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
    {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
    {{0.0, 10.0, 100.0, 200.0}, {0.0, 10.0, 100.0, 200.0}, {1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0}},
    {{0.0, 1.0, 0.0, 1.0}, {0.0, 0.0, 1.0, 2.0}, {2.0, 10.0, 0.0, 5.}, {2, -2.56637, 0, -1.28318531}}

};

void test_simulation()
{

    for (int i = 0; i < sizeof(test_cases) / sizeof(test_case); i++)
    {
        std::vector<double> x = test_cases[i].x;
        std::vector<double> y = test_cases[i].y;
        std::vector<double> theta = test_cases[i].theta;
        std::vector<double> correct_result = test_cases[i].correct_result;

        std::cout << "Testing simulation function: " << i << std::endl;
        std::vector<double> res = simulation(x, y, theta, R);

        assert(res.size() == correct_result.size());
        for (int i = 0; i < res.size(); i++)
        {
            // print type of res[i] and correct_result[i]
            double diff = res[i] - correct_result[i];
            if (diff > threshold || diff < -threshold)
            {
                std::cout << "res[i]: " << res[i] << " correct_result[i]: " << correct_result[i] << std::endl;
                std::cout << "diff: " << diff << std::endl;
                assert(res[i] == correct_result[i]);
            }
        }
    }
}

void test_simulation_openmp()
{

    for (int i = 0; i < sizeof(test_cases) / sizeof(test_case); i++)
    {
        std::vector<double> x = test_cases[i].x;
        std::vector<double> y = test_cases[i].y;
        std::vector<double> theta = test_cases[i].theta;
        std::vector<double> correct_result = test_cases[i].correct_result;

        std::cout << "Testing simulation_openmp function: " << i << std::endl;
        std::vector<double> res = simulation_openmp(x, y, theta, R);

        assert(res.size() == correct_result.size());
        for (int i = 0; i < res.size(); i++)
        {
            // print type of res[i] and correct_result[i]
            double diff = res[i] - correct_result[i];
            if (diff > threshold || diff < -threshold)
            {
                std::cout << "res[i]: " << res[i] << " correct_result[i]: " << correct_result[i] << std::endl;
                std::cout << "diff: " << diff << std::endl;
                assert(res[i] == correct_result[i]);
            }
        }
    }
}

void test_simulation_mpi()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get each process rank

    for (int i = 0; i < sizeof(test_cases) / sizeof(test_case); i++)
    {
        std::vector<double> x = test_cases[i].x;
        std::vector<double> y = test_cases[i].y;
        std::vector<double> theta = test_cases[i].theta;
        std::vector<double> correct_result = test_cases[i].correct_result;
        if (rank == 0)
            std::cout << "Testing simulation_mpi function: " << i << std::endl;

        std::vector<double> res = simulation_mpi(x, y, theta, R);

        if (rank == 0)
        { // Only assert the result on the root process
            assert(res.size() == correct_result.size());
            for (int i = 0; i < res.size(); i++)
            {
                // print type of res[i] and correct_result[i]
                double diff = res[i] - correct_result[i];
                if (diff > threshold || diff < -threshold)
                {
                    std::cout << "res[i]: " << res[i] << " correct_result[i]: " << correct_result[i] << std::endl;
                    std::cout << "diff: " << diff << std::endl;
                    assert(res[i] == correct_result[i]);
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int provided, rank;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get each process rank

    if (rank == 0)
    {
        std::cout << "Running tests..." << std::endl;
    }

    if (rank == 0)
    {
        test_simulation();
        test_simulation_openmp();
    }

    test_simulation_mpi();

    if (rank == 0)
        std::cout << "All tests passed!" << std::endl;

    MPI_Finalize();
}
