#include "namespace.h"
#include "util.h"

using namespace std;

long double computeEuclideanDist(std::vector<long double>& point_i, std::vector<long double>& point_j)
{
	long double dist = 0.0;
	for (int dd=0 ; dd < point_i.size() ; ++dd)
	{
		dist += (point_i[dd]-point_j[dd])*(point_i[dd]-point_j[dd]);
	}
	return sqrt(dist);
}


bool sortcol(const vector<long double>& v1, const vector<long double>& v2) 
{
	return v1[0] < v2[0];
}


int printVec(vector<long double> &a)
{
	// cout << "Printing point -- " << endl;
	for(int i = 0; i < a.size(); ++i)
	{
		cout << a[i] << "\t";
	}

	cout << endl;

	return 0;
}


int printVecVec(vector <vector<long double>> &p)
{
	for(int i = 0; i < p.size(); i++)
	{
		cout << "point Id: " << i << " ---- ";
		printVec(p[i]);
	}

	return 0;
}

vector<string> getPathAndFileName(string fileName)
{
	vector<string> pathAndFileName;

	std::size_t found = fileName.find_last_of("/");
	pathAndFileName.push_back(fileName.substr(0,found));
	pathAndFileName.push_back(fileName.substr(found+1));

	return pathAndFileName;
}