# Parallelizing the Vicsek Model

## Getting Started
To compile and run the main simulation:

```bash
make run
```
This command will build the main simulation binary and execute it.
Build the Project

To compile the standard version of the simulation:

```bash

make birds
```
This will compile the main simulation program without OpenMP support.

To compile the OpenMP-enabled version of the simulation:

```bash

make openMP
```
This will compile the main simulation program with OpenMP support.


```bash
make mpi
```
This will compile the main simulation program with MPI support.

To compile and run the tests:

```bash
make test
```
This will compile the test program and execute it.
Benchmark the Project

To compile the benchmarking program:

```bash
make bench
```
This will produce a binary named bench.out located in the bin folder.

To remove all compiled binary and object files:
```bash
make clean
```
This will clean up the bin directory by removing all object files and executables.


## View Documentation
To view the documentation, open the file ``files.html`` in the directory named ``html``.
