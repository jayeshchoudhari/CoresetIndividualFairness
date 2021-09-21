#ifndef UTIL_H
#define UTIL_H

long double computeEuclideanDist(std::vector<long double>& point_i, std::vector<long double>& point_j);
bool sortcol(const std::vector<long double>& v1, const std::vector<long double>& v2);
int printVec(std::vector<long double> &a);
int printVecVec(std::vector <std::vector<long double>> &p);
std::vector<std::string> getPathAndFileName(std::string fileName);

#endif // UTIL_H