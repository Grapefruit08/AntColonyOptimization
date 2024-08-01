#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>

#include <chrono>

using namespace std;

class AntColonyOptimization_sequential {
private:
	int numNodes;
	int numAnts;
	double alpha;
	double beta;
	double evaporation;
	double Q;

	double bestLength = DBL_MAX;

	vector<int> ConstructTour(int start = -1);
	int ChooseNextCity(vector<int>& tour, vector<bool>& visited);
	void UpdatePheromone(const vector<vector<int>>& ant_tours, const vector<double>& tour_lengths);
	void loadDistanceMatrix(const std::string& filename);

public:
	vector<vector<double>> distance_matrix;
	vector<vector<double>> pheromone_matrix;

	AntColonyOptimization_sequential(int num_cities, int num_ants, double alpha, double beta, double evaporation, double Q) :
		numNodes(num_cities), numAnts(num_ants), alpha(alpha), beta(beta), evaporation(evaporation), Q(Q) {

		pheromone_matrix.assign(num_cities, vector<double>(num_cities, 1.0));
		distance_matrix.assign(num_cities, vector<double>(num_cities));

		// initialize random lengths between cities, symmetrically
		for (int i = 0; i < num_cities; ++i) {
			for (int j = i + 1; j < num_cities; ++j) {
				double dist = rand() % 1000 + 1; // random distance between 1 and 100
				distance_matrix[i][j] = distance_matrix[j][i] = dist;
			}
		}
	};

	AntColonyOptimization_sequential(const std::string& filename, int num_ants, double alpha, double beta, double evaporation, double Q)
		: numAnts(num_ants), alpha(alpha), beta(beta), evaporation(evaporation), Q(Q) {

		// Load the distance matrix from the file
		loadDistanceMatrix(filename);

		// Initialize the pheromone matrix
		pheromone_matrix.assign(numNodes, std::vector<double>(numNodes, 1.0));
	}

	vector<int> RunACO(int max_iterations);

	double CalcLengthTour(const vector<int>& tour);

	void PrintTour(const vector<int>& tour, string msg);

	void SetDistMatrix(const vector<vector<double>>& new_distance_matrix);

	void WriteMatrixToFile(const vector<vector<double>>& matrix, const string& filename);
};
