# AntColonyOptimization

## Overview
This project implements an Ant Colony Optimization (ACO) algorithm to solve the Traveling Salesman Problem (TSP) using both single-threaded and multi-threaded (MPI) approaches. The ACO algorithm is a probabilistic technique inspired by the foraging behavior of ants to find approximate solutions to combinatorial optimization problems.

## Files

 - **AntColonyOptimization.h**:  Header file for the single-threaded ACO implementation.
 - **AntColonyOptimization.cpp**: Contains the single-threaded implementation of the Ant Colony Optimization algorithm.
 - **AntColonyOptimization_MPI.h**: Header file for the MPI implementation.
 - **AntColonyOptimization_MPI.cpp**: Contains the multi-threaded (MPI) implementation of the Ant Colony Optimization algorithm.
 - **Main.cpp**: Contains the main function to run the ACO algorithm.
 - **dist_*N*.txt**: Contains the distance matrix for a *N*-city TSP problem.

## Requirements

- C++11 or higher
- MPI library (for the multi-threaded version)

## Usage

### Input

- The program reads the distance matrix from e.g. `dist_10.txt`. Ensure this file is in the same directory as the executable.
- The format of `dist_10.txt` is a matrix where the entry in the i-th row and j-th column represents the distance between city i and city j.

### Output

- The program will output the best path found and its total distance.

## Code Structure

### AntColonyOptimization.h
This header file declares the `AntColonyOptimization_sequential` class, which implements the single-threaded ACO algorithm. It includes:

- **ConstructTour**: Constructs a tour for an ant.
- **ChooseNextCity**: Chooses the next city for an ant to visit.
- **UpdatePheromone**: Updates the pheromone levels based on the tours of all ants.
- **loadDistanceMatrix**: Loads the distance matrix from a file.
- **RunACO**: Runs the ACO algorithm for a specified number of iterations.
- **CalcLengthTour**: Calculates the total length of a given tour.
- **PrintTour**: Prints the tour.
- **SetDistMatrix**: Sets the distance matrix.
- **WriteMatrixToFile**: Writes a matrix to a file.

### AntColonyOptimization.cpp
This file contains the implementation of the methods declared in AntColonyOptimization.h for solving the TSP in a single-threaded manner.

### AntColonyOptimization_MPI.h and AntColonyOptimization_MPI.cpp
These files contain the implementation of the ACO algorithm for solving the TSP using MPI for parallel processing. The header file declares the functions and the MPI-specific data structures used in the MPI implementation.

### Main.cpp
This file contains the main function, which initializes the parameters for the ACO algorithm, reads the distance matrix from .txt file, and starts the optimization process by calling the appropriate functions.

### dist_*N*.txt
This file contains the distance matrix for a *N*-city TSP problem. The distances are stored in a plain text format, with rows and columns representing the distances between different cities.



