#include "AntColonyOptimization_MPI.h"

void AntColonyOptimization::loadDistanceMatrix(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> temp_matrix;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double distance;
        std::vector<double> row;
        while (ss >> distance) {
            row.push_back(distance);
        }
        temp_matrix.push_back(row);
    }

    file.close();

    // Set the number of nodes
    numNodes = temp_matrix.size();

    // Assign the loaded matrix to the distance matrix
    distance_matrix = temp_matrix;
};

void AntColonyOptimization::RunAco(int numIterations)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int iter = 0; iter < numIterations; iter++) {

        vector<int> oneProcessBestTour;
        double oneProcessBestLength = numeric_limits<double>::max();

        vector<vector<int>> oneProcessTours;
        vector<double> oneProcessTourLengths;

        // distribute ants to the proccesses
        for (int ant = rank; ant < numAnts; ant += size) {
            vector<int> tour = ConstructTour(-1);
            double tourLength = CalcLengthTour(tour);

            oneProcessTours.push_back(tour);
            oneProcessTourLengths.push_back(tourLength);

            if (tourLength < oneProcessBestLength) {
                oneProcessBestTour = tour;
                oneProcessBestLength = tourLength;
            }
        }

        // Find the best length for the iter-th iteration
        double iterationBestTourLength;
        MPI_Allreduce(&oneProcessBestLength, &iterationBestTourLength, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        // if iter-th bestTourLength smaller than all time bestTourLength -> comminication
        // otherwise skip to the pheromone update
        if (iterationBestTourLength < bestTourLength) {
            // find the proccesse with the shortest path
            int best_rank;
            if (oneProcessBestLength == iterationBestTourLength) {
                best_rank = rank;
            }
            else {
                best_rank = size; // if the proccesse doesnt find the shortest path set his best_rank to size
                // bellow MIN Reduction will be done over best_rank
            }

            // chose the proccesse with the lowest best_rank in case more 
            // than one proccesse found the shoertest path of the iteration
            int global_best_rank;
            MPI_Allreduce(&best_rank, &global_best_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

            // update bestTour in all proccesses from the global_best_rank proccesse
            if (rank == global_best_rank) {
                bestTour = oneProcessBestTour;
            }
            MPI_Bcast(bestTour.data(), numNodes, MPI_INT, global_best_rank, MPI_COMM_WORLD);

            // set new bestTourLength
            bestTourLength = iterationBestTourLength;
        }

        std::vector<std::vector<double>> oneProcessAddition =
            ComputePheromoneAddition(oneProcessTours, oneProcessTourLengths);

        // pheromone update, in rank 0 evaporate, set the last pheromone_matrix as its addition
        if (rank == 0) {
            EvaporatePheromone();
            //print_all_data(&rank);
            oneProcessAddition = AddPheromone(oneProcessAddition);
        }

        //print_all_data(&rank);


        // set new pheromone matrix based on additions of all
        for (int i = 0; i < numNodes; i++) {
            MPI_Allreduce(oneProcessAddition[i].data(), pheromone_matrix[i].data(), numNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
    }
}

void AntColonyOptimization::printResults() {
    std::cout << "Best tour: ";
    for (int i = 0; i < numNodes; i++) {
        std::cout << bestTour[i] << " ";
    }
    std::cout << "\nTour length: " << bestTourLength << std::endl;
}

void AntColonyOptimization::printTour(std::vector<int>& tour, int tourLength) {
    for (int i = 0; i < numNodes; i++) {
        std::cout << tour[i] << " ";
    }
    std::cout << "\nTotal length: " << tourLength << std::endl;
}

double AntColonyOptimization::CalcLengthTour(const std::vector<int>& tour) {
    double total_length = 0.0;
    for (int i = 1; i < numNodes; i++) {
        int prev_city = tour[i - 1];
        total_length += distance_matrix[prev_city][tour[i]];
    }

    // return to the starting city
    total_length += distance_matrix[tour[numNodes - 1]][tour[0]];

    return total_length;
}

std::vector<int> AntColonyOptimization::ConstructTour(int start) {
    vector<int> tour;
    vector<bool> visited(numNodes, false);

    // choose the starting city randomly
    int start_city = start;
    if (start == -1) {
        // better generator than rand() // random number from <0,1>
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, numNodes - 1);

        //start_city = rand() % numNodes;
        start_city = dis(gen);
    }

    // add start_city as the first of visited cities 
    tour.push_back(start_city);
    visited[start_city] = true;

    // add all cities to the tour
    for (int i = 1; i < numNodes; i++) {
        int next_city = ChooseNextCity(tour, visited);
        tour.push_back(next_city);
        visited[next_city] = true;
    }

    return tour;
};

int AntColonyOptimization::ChooseNextCity(std::vector<int>& tour, std::vector<bool>& visited) {
    int current_city = tour.back();

    double total_probability = 0.0;
    vector<double> probabilities(numNodes, 0.0);

    // calculate probability for choosing i-th city
    for (int i = 0; i < numNodes; i++) {
        if (visited[i] == false) {
            double pheromone_factor = pow(pheromone_matrix[current_city][i], alpha);
            double heuristic_factor = pow(1.0 / distance_matrix[current_city][i], beta);
            double probability = pheromone_factor * heuristic_factor;

            probabilities[i] = probability;
            total_probability += probability;
        }
    }

    // better generator than rand() // random number from <0,1>
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    //double random_number = (double)rand() / RAND_MAX; // random number from <0,1>
    double random_number = dis(gen); // random number from <0,1>
    double cumulated_probab = 0.0;
    for (int i = 0; i < numNodes; i++) {
        if (visited[i] == false) {
            cumulated_probab += probabilities[i] / total_probability;
            if (cumulated_probab >= random_number) {
                return i;
            }
        }
    }

    return -1; // this means problem and function should never come reach here
}

void AntColonyOptimization::EvaporatePheromone() {
    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            pheromone_matrix[i][j] *= (1.0 - evaporation);
            pheromone_matrix[j][i] = pheromone_matrix[i][j]; // Pheromone matrix is symetric
        }
    }
}

std::vector<std::vector<double>> AntColonyOptimization::AddPheromone(const std::vector<std::vector<double>>& pheromoneAddition) {
    for (int i = 0; i < numNodes; i++) {
        for (int j = i; j < numNodes; j++) {
            pheromone_matrix[i][j] += pheromoneAddition[i][j];
            pheromone_matrix[j][i] = pheromone_matrix[i][j];  // Pheromone matrix is symetric
        }
    }
    return pheromone_matrix;
}

std::vector<std::vector<double>> AntColonyOptimization::ComputePheromoneAddition
(std::vector<std::vector<int>>& oneProcessTours, std::vector<double>& oneProcessTourLengths)
{
    std::vector<std::vector<double>> pheromoneAdd(numNodes, vector<double>(numNodes, 0.0));

    for (int ant = 0; ant < oneProcessTourLengths.size(); ant++) {
        double the_ants_pheromone = Q / oneProcessTourLengths[ant];
        for (int i = 1; i < numNodes; i++) {
            int city_1 = oneProcessTours[ant][i - 1];
            int city_2 = oneProcessTours[ant][i];

            pheromoneAdd[city_1][city_2] += the_ants_pheromone;
            pheromoneAdd[city_2][city_1] += the_ants_pheromone;
        }

        // add pheromon also to the edge connecting last node with starting one
        int city_1 = oneProcessTours[ant][numNodes - 1];
        int city_2 = oneProcessTours[ant][0];

        pheromoneAdd[city_1][city_2] += the_ants_pheromone;
        pheromoneAdd[city_2][city_1] += the_ants_pheromone;
    }
    return pheromoneAdd;
};

// just testing if sending is ok
void AntColonyOptimization::print_all_data(int* rank) {
    std::cout << "Rank: " << *rank << ": numNodes = " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++) {
        std::cout << *rank << ": ";
        for (int j = 0; j < numNodes; j++) {
            std::cout << pheromone_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Rank: " << *rank << std::endl;
}
