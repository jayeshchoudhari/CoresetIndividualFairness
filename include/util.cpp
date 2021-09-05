#include "namespace.h"
#include "util.h"

using namespace std;

long double computeEuclideanDist(std::vector<long double> point_i, std::vector<long double> point_j)
{
	long double dist = 0.0;
	for (int dd=0 ; dd < point_i.size() ; ++dd)
	{
		dist += (point_i[dd]-point_j[dd])*(point_i[dd]-point_j[dd]);
	}
	return sqrt(dist);
}