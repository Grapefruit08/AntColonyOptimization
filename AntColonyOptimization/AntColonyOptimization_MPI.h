#pragma once

#include <mpi.h>
#include <vector>
#include <iostream>
#include <random>
#include <sstream>
#include <fstream>

using namespace std;

class AntColonyOptimization {
private:
    std::vector<std::vector<double>> distance_matrix;
    std::vector<std::vector<double>> pheromone_matrix;
    std::vector<int> bestTour;
    double bestTourLength;

    int numNodes;
    int numAnts;
    double alpha;
    double beta;
    double evaporation;
    double Q;

    double CalcLengthTour(const std::vector<int>& tour);

    std::vector<int> ConstructTour(int start);

    int ChooseNextCity(std::vector<int>& tour, std::vector<bool>& visited);
    void EvaporatePheromone();

    std::vector<std::vector<double>> AddPheromone(const std::vector<std::vector<double>>& pheromoneAddition);

    std::vector<std::vector<double>>ComputePheromoneAddition
    (std::vector<std::vector<int>>& oneProcessTours, std::vector<double>& oneProcessTourLengths);

    void loadDistanceMatrix(const std::string& filename);

public:
    void RunAco(int numIterations);
    void printResults();

    void printTour(std::vector<int>& tour, int tourLength);

    void print_all_data(int* rank); // jenom testovaci vypis

    AntColonyOptimization(const std::string& filename, int num_ants, double alpha, double beta, double evaporation, double Q)
        : numAnts(num_ants), alpha(alpha), beta(beta), evaporation(evaporation), Q(Q) {

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        bestTourLength = std::numeric_limits<double>::max();

        if (rank == 0) {
            // Load the distance matrix from the file
            loadDistanceMatrix(filename);

            // Send numNodes and distance matrix to other processes
            for (int proc = 1; proc < size; ++proc) {
                MPI_Send(&numNodes, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
                for (int j = 0; j < numNodes; ++j) {
                    MPI_Send(distance_matrix[j].data(), numNodes, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
                }
            }
        }
        else {
            // Receive numNodes from the root process
            MPI_Recv(&numNodes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Resize the distance matrix
            distance_matrix.resize(numNodes, std::vector<double>(numNodes));

            // Receive the distance matrix from the root process
            for (int i = 0; i < numNodes; ++i) {
                MPI_Recv(distance_matrix[i].data(), numNodes, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // Initialize the pheromone matrix and bestTour
        pheromone_matrix.assign(numNodes, std::vector<double>(numNodes, 1.0));
        bestTour.assign(numNodes, (0));
    }
};
