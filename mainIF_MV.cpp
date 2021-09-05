#include "./include/namespace.h"
#include "./include/dataIO.h"

const int maxn = 11000;
const int maxd = 30;
const long double eps = 1e-6;

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 5) 
	{
		std::cerr << "ERROR -- Requires 4 parameters: 1) Input File name; 2) Number of points (n); 3) Dimension of data (d); and 4) number of centers (k).\n";
        exit(1);
	}

	string inputFileName = argv[1];
	int n = atoi(argv[2]);
    int d = atoi(argv[3]);
    int k = atoi(argv[4]);

    DataIO DIO(n, d, k, inputFileName, maxn, maxd, eps);

    DIO.computeFairRadius(maxn);

    DIO.computeJungEtAlAlphaAndCost();

    DIO.computeMVIndFair();

    DIO.generate_output();

	return 0;
}