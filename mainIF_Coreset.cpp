#include "./include/namespace.h"
#include "./include/dataIO.h"
#include <chrono>

const int maxn = 10000;
const int maxd = 10;
const long double eps = 1e-6;

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
        // std::cerr << "ERROR -- Requires 6 parameters: 1) Input File name; 2) Number of points (n); 3) Dimension of data (d); 4) number of centers (k); 5)Full data filename; and 6) Fair Radius from data (file) from MV ...\n";
		std::cerr << "ERROR -- Requires 4 parameters: 1) Input File name; 2) number of centers (k); 3)Full data filename; and 4) Fair Radius from data (file) from MV ...\n";
        exit(1);
	}

	string inputFileName = argv[1];
	// int n = atoi(argv[2]);
    // int d = atoi(argv[3]);
    int k = atoi(argv[2]);
    string fullDataFileName = argv[3];
    string fairRadiusFileName = argv[4];

    // DataIO DIO(n, d, k, inputFileName, maxn, maxd, eps);
    DataIO DIO(k, inputFileName, maxn, maxd, eps);
    DIO.readInputWeighted(inputFileName);

    //time start
    auto startTime = std::chrono::system_clock::now();

    cout << "Computing weighted Fair radius...\n";
    DIO.computeWeightedFairRadius(maxn);
    cout << "Got weighted Fair radius...n";

    DIO.computeJungEtAlAlphaAndCostWeighted();

    DIO.computeMVIndFairWeighted();

    auto endTime = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    double totalComputeTime = elapsed.count();
    //time end

    DIO.computeCostOnFullDataUsingCoresetCenters(fullDataFileName);

    DIO.computeFairnessOnFullDataUsingCoresetCenters(fairRadiusFileName);

    DIO.generate_output(inputFileName, totalComputeTime, 1);

	return 0;
}