#include "./include/namespace.h"
#include "./include/dataIO.h"
#include <chrono>

const int maxn = 10000;
const int maxd = 10;
const long double eps = 1e-6;

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 3) 
	{
        // std::cerr << "ERROR -- Requires 4 parameters: 1) Input File name; 2) Number of points (n); 3) Dimension of data (d); and 4) number of centers (k).\n";
		std::cerr << "ERROR -- Requires 2 parameters: 1) Input File name; 2) number of centers (k).\n";
        exit(1);
	}

	string inputFileName = argv[1];
	// int n = atoi(argv[2]);
    // int d = atoi(argv[3]);
    int k = atoi(argv[2]);

    // DataIO DIO(n, d, k, inputFileName, maxn, maxd, eps);
    DataIO DIO(k, inputFileName, maxn, maxd, eps);
    DIO.readInput(inputFileName);

    //time start
    auto startTime = std::chrono::system_clock::now();
    DIO.computeFairRadius(maxn, inputFileName);

    DIO.computeJungEtAlAlphaAndCost();

    DIO.computeMVIndFair();
    
    auto endTime = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    double totalComputeTime = elapsed.count();

    DIO.generate_output(inputFileName, totalComputeTime, 0);

	return 0;
}