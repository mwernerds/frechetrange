// Entrypoint for the program, loads queryset and dataset and invokes main algorithm
// Also determines number of worker threads used by the algorithm
#include "tue.hpp"
#include <stdio.h>
#include <string>
#include <thread>

using namespace tue_details;

int main(int argc, char *argv[])
{
	std::string datasetFilename = "dataset.txt";
	std::string querysetFilename = "queries.txt";
	if (argc == 3) {
		datasetFilename = std::string(argv[1]);
		querysetFilename = std::string(argv[2]);
	}

	long timeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	
	std::cout << "Dataset: " << datasetFilename << " Queryset: " << querysetFilename << "\n";

	BoundingBox *box = new BoundingBox();

	AlgoData a;
	a.queries = a.fio.parseQueryFile(querysetFilename);
	std::cout << "Loaded queries\n";
	a.trajectoryNames = a.fio.parseDatasetFile(datasetFilename);
	a.numTrajectories = a.trajectoryNames->size();
	std::cout << "Loaded trajectories\n";
	a.numWorkers = std::thread::hardware_concurrency();// worker threads == number of logical cores
	a.boundingBox = box;

	#if !USE_MULTITHREAD
		a.numWorkers = 1;
	#endif


	std::cout << "Loaded dataset and query files\n";

	std::cout << "Num workers: " << a.numWorkers << "\n";

	runAlgorithm(&a);

	timeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1) - timeMS;

	std::cout << "Total time: " << timeMS << "\n";


	std::cout << "Done\n";



	return 0;
}
