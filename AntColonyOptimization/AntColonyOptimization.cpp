#include "AntColonyOptimization.h"

vector<int> AntColonyOptimization_sequential::ConstructTour(int start) {
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

int AntColonyOptimization_sequential::ChooseNextCity(vector<int>& tour, vector<bool>& visited) {
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
};

double AntColonyOptimization_sequential::CalcLengthTour(const vector<int>& tour) {
	double total_length = 0.0;
	for (int i = 1; i < numNodes; i++) {
		int prev_city = tour[i - 1];
		total_length += distance_matrix[prev_city][tour[i]];
	}

	// return to the starting city
	total_length += distance_matrix[tour[numNodes - 1]][tour[0]];

	return total_length;
};

void AntColonyOptimization_sequential::UpdatePheromone(const vector<vector<int>>& ant_tours, const vector<double>& tour_lengths) {
	// evaporate phoromone
	for (int i = 0; i < numNodes; ++i) {
		for (int j = i + 1; j < numNodes; ++j) {
			pheromone_matrix[i][j] *= (1.0 - evaporation);
			pheromone_matrix[j][i] = pheromone_matrix[i][j]; // Pheromone matrix is symetric
		}
	}

	for (int ant = 0; ant < numAnts; ant++) {
		double the_ants_pheromone = Q / tour_lengths[ant];
		for (int i = 1; i < numNodes; i++) {
			int city_1 = ant_tours[ant][i - 1];
			int city_2 = ant_tours[ant][i];

			pheromone_matrix[city_1][city_2] += the_ants_pheromone;
			pheromone_matrix[city_2][city_1] += the_ants_pheromone;
		}

		// oferemovat i cestu zpet do start_city
		int city_1 = ant_tours[ant][numNodes - 1];
		int city_2 = ant_tours[ant][0];

		pheromone_matrix[city_1][city_2] += the_ants_pheromone;
		pheromone_matrix[city_2][city_1] += the_ants_pheromone;
	}
}
void AntColonyOptimization_sequential::loadDistanceMatrix(const std::string& filename) {
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

vector<int> AntColonyOptimization_sequential::RunACO(int max_iterations) {
	vector<int> best_tour;
	//double best_length = INF;

	for (int i = 0; i < max_iterations; i++) {
		vector<vector<int>> ant_tours;
		vector<double> tour_lengths;

		// release the ants
		for (int ant = 0; ant < numAnts; ant++) {
			vector<int> ant_tour;
			if (numAnts == numNodes) {
				ant_tour = ConstructTour(ant); // for smaller graphs faster konvergation
			}
			else {
				ant_tour = ConstructTour();
			}

			double tour_length = CalcLengthTour(ant_tour);

			//PrintTour(ant_tour, tour_length, "iteration: " + std::to_string(i));

			ant_tours.push_back(ant_tour);
			tour_lengths.push_back(tour_length);

			// save the best tour
			if (tour_lengths[ant] < bestLength) {
				best_tour = ant_tour;
				bestLength = tour_length;
			}
		}

		UpdatePheromone(ant_tours, tour_lengths);
	}

	return best_tour;
}
void AntColonyOptimization_sequential::PrintTour(const vector<int>& tour, string msg) {
	if (msg != "") {
	std:cout << msg << std::endl;
	}
	std::cout << "Start -> ";
	for (int i = 0; i < numNodes; i++) {
		std::cout << tour[i] << " ";
	}
	std::cout << " -> End, Length: " << bestLength << std::endl;
};

void AntColonyOptimization_sequential::SetDistMatrix(const vector<vector<double>>& new_distance_matrix) {
	if (new_distance_matrix.size() != numNodes || new_distance_matrix[0].size() != numNodes) {
		cout << "Error: Given matrix size does not match the number of cities." << endl;
		return;
	}

	for (int i = 0; i < numNodes; ++i) {
		for (int j = 0; j < numNodes; ++j) {
			distance_matrix[i][j] = new_distance_matrix[i][j];
		}
	}
}

void AntColonyOptimization_sequential::WriteMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
	ofstream outFile(filename);

	if (!outFile.is_open()) {
		cout << "Error opening file " << filename << endl;
		return;
	}

	int width = 8;

	// Write the matrix
	for (const auto& row : matrix) {
		for (double matrix_row : row) {
			outFile << setw(width) << matrix_row << "\t"; // tab-separated
		}
		outFile << endl;
	}

	outFile.close();
}
