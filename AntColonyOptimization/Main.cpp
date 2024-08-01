#include "AntColonyOptimization.h"
#include "AntColonyOptimization_MPI.h"

int main(int argc, char** argv) {

    srand(123);

    double alpha = 1.0;
    double beta = 5.0;
    double evaporation = 0.1;
    double Q = 100;

    int numAnts = 20;
    int numIterations = 100;

    string matrix_file = "dist_50.txt";

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Sequential ACO (ACO_1) - Only on rank 0
    //-----------------------------------------------------------------------//
    if (rank == 0) {
        std::cout << "Start ACO_sequential" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        AntColonyOptimization_sequential aco(matrix_file, numAnts, alpha, beta, evaporation, Q);
        std::vector<int> best_tour = aco.RunACO(numIterations);
        aco.PrintTour(best_tour, "Best Tour: ");

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Time taken by ACO_sequential: " << duration.count() << " seconds" << std::endl;
    }

    // Synchronize processes before starting MPI ACO
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI ACO (ACO_MPI)
    //-----------------------------------------------------------------------//
    {
        double start, end, duration;
        if (rank == 0) {
            std::cout << "________________________________________________________________" << std::endl;
            std::cout << "Start ACO_MPI" << std::endl;
            start = MPI_Wtime();
        }

        // Ant Colony Optimization setup
        AntColonyOptimization aco(matrix_file, numAnts, alpha, beta, evaporation, Q);
        aco.RunAco(numIterations);

        // Print final result of ACO
        if (rank == 0) {
            aco.printResults();
        }

        // Synchronize processes before measuring the end time
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            end = MPI_Wtime();
            duration = end - start;
            std::cout << "Time taken by ACO_MPI: " << duration << " seconds" << std::endl;
        }
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
