#include "namespace.h"
#include "dataIO.h"
#include "util.h"

using namespace std;

// DataIO :: DataIO(int numPoints, int dim, int kcenters, string inputFileName, int maxn, int maxd, long double epsVal)
DataIO :: DataIO(int kcenters, string inputFileName, int maxn, int maxd, long double epsVal)
{
	// n = numPoints;
	// d = dim;
	k = kcenters;
	eps = epsVal;

	// cout << n << " " << d  << " " << k << " " << eps << endl;
	cout << k << " " << eps << endl;

	// point.resize(maxn, vector<long double>(maxd));
	// readInput(inputFileName);
}


//Reading the input points. Each point comes in a separate line and the coordinates are separated by space.
int DataIO :: readInput(string fileName)
{
	stringstream ss;
	string line;

	ifstream fin(fileName);
	// space separated  values...
	int i = 0;
	while(getline(fin, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		vector<long double> tempVec;
		long double tempval;

		while(ss >> tempval)
		{
			tempVec.push_back(tempval);
		}
		
		// cout << "push_back point ... \n";
		point.push_back(tempVec);

		corWeights.push_back(1.0);
		// cout << "point ... \n";

		i += 1;
	}

	fin.close();

	n = point.size();
	d = point[0].size();

	cout << "Going for initialization ... \n";
	initializeDataStructure(n);

	return 0;
}

/*
//Reading the input points. Each point comes in a separate line and the coordinates are separated by space.
int DataIO :: readInput(string fileName)
{
	cout << "Reading input -- " << endl;
	cout << n << " " << d  << " " << k << " " << eps << endl;

	long double tempVar;
	// ifstream fin("input.csv");
	ifstream fin(fileName);
	// space separated  values...
	for(int i = 0; i < n; i++)
	{
		vector<long double> tempvec;
		long double tempval;
		for(int j = 0; j < d; j++)
		{
			fin >> tempval;
			tempvec.push_back(tempval);
		}
		point.push_back(tempvec);

		// if there is any other with weights.. 
		// last index of each line is weight...
		corWeights[i] = 1.0;
	}

	fin.close();
	// printVecVec(point);
	return 0;
}
*/


int DataIO :: readInputWeighted(string fileName)
{
	stringstream ss;
	string line;

	ifstream fin(fileName);
	// space separated  values...
	int i = 0;
	long double weightVal;
	while(getline(fin, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		vector<long double> tempVec;
		long double tempval;

		while(ss >> tempval)
		{
			tempVec.push_back(tempval);
		}
		// if there is any other with weights...
		// last index of each line is weight...
		weightVal = tempVec.back();
		corWeights.push_back(weightVal);

		// popping back the last weight element 
		tempVec.pop_back();
		point.push_back(tempVec);

		i += 1;
	}

	fin.close();
	
	n = point.size();
	d = point[0].size();
	
	initializeDataStructure(n);

	return 0;
}


/*
// Reading the input points. Each point comes in a separate line and the coordinates are separated by comma.
// last index is the weight of the point
int DataIO :: readInputWeighted(string fileName)
{
	cout << "Reading input -- " << endl;
	cout << n << " " << d  << " " << k << " " << eps << endl;

	ifstream fin(fileName);
	// space separated  values...
	for(int i = 0; i < n; i++)
	{
		long double weightVal;
		vector<long double> tempvec;
		long double tempval;
		for(int j = 0; j < d; j++)
		{
			fin >> tempval;
			tempvec.push_back(tempval);
		}
		point.push_back(tempvec);
		fin >> weightVal;

		// if there is any other with weights...
		// last index of each line is weight...
		corWeights[i] = weightVal;
	}

	fin.close();
	// printVecVec(point);
	return 0;
}
*/

int DataIO :: initializeDataStructure(int numPoints)
{
	cout << "Starting to initialize\n";
	r.resize(numPoints);
	// ogWeights.resize(maxn, 0.0);
	// corWeights.resize(numPoints, 0.0);

	covered.resize(numPoints, 0);	//shows if a point is covered by current centers
	center.resize(numPoints, -1);	//shows the corresponding center for each point -- default set to -1

	whichBall.resize(numPoints, -1);	//keeps for each point in which critical ball they exist or -1 otherwise
	iscenter.resize(numPoints, 0);	//keeps for each point if it is in the current set of centers

	cout << "Initialized Data-structures -- \n";

	return 0;
}


long double DataIO :: compute_r(int i, int maxn)
{
	//computes the fair radius for each point
	long double temp[n];

	for (int j=0 ; j<n ; ++j)
	{
		// temp[j] = computeEuclideanDist(point[i], point[j]);
		temp[j] = compute_dist(i, j);
		// cout << i << " " << temp[j] << endl;
	}

	sort(temp, temp+n);

	// cout << "i = " << i << " rVal : " << temp[(n-1)/k] << " (n-1)/k = " << (n-1)/k << endl;

	// have to consider weights... 
	// take sum of first l points that weigh to n-1/k -- change
	return temp[(n-1)/k];
}


long double DataIO :: compute_r_weighted(int i, int maxn, long double totalWeight)
{
	//computes the fair radius for each point
	long double temp[n];
	vector< vector<long double> > tempRadWeights;
	// tempRadWeights.resize(maxn, vector<long double>(2));

	for (int j=0 ; j<n ; ++j)
	{
		// temp[j] = computeEuclideanDist(point[i], point[j]);
		temp[j] = compute_dist(i, j);
		vector<long double>  tempVec;

		tempVec.push_back(temp[j]);
		tempVec.push_back(corWeights[j]);

		tempRadWeights.push_back(tempVec);
	}

	// sort(tempRadWeights.begin(), tempRadWeights.end());
	sort(tempRadWeights.begin(), tempRadWeights.end(), sortcol);

	// printVecVec(tempRadWeights);

	long double localWeight = 0.0;
	// This is over expectation -- n is not present...
	// Change this to sum of weights...
	// long double maxWeight = ((n-1)*1.0)/k;
	long double maxWeight = totalWeight/k;
	long double rVal = tempRadWeights[(n-1)/k][0];
	// cout << "max weight : " << maxWeight << " totalWeight = " << totalWeight << " initial rval = " << rVal << " (n-1)/k = " << (n-1)/k <<  endl;

	for(int q = 0; q < tempRadWeights.size(); q++)
	{
		if(localWeight + tempRadWeights[q][1] <= maxWeight)
		{
			localWeight += tempRadWeights[q][1];
		}
		else if(q == 0)
		{
			return tempRadWeights[0][0];
		}
		else
		{
			rVal = tempRadWeights[q-1][0];
			break;
		}
	}

	// cout << "i = " << i << " rVal : " << rVal << endl;
	// have to consider weights... 
	// take sum of first l points that weigh to n-1/k -- change
	// return temp[(n-1)/k];
	return rVal;
}

int DataIO :: writeFairRadiusToFile(string fileName)
{
	string strk = std::to_string(k);

	vector<string> pathFileName = getPathAndFileName(fileName);

	string fairOutFileName = "fairRadiusFiles/fairRadius_" + strk + "_" + pathFileName[1];
	fairOutFileName = pathFileName[0] + "/" + fairOutFileName;

	cout << "writing to fairOutFileName\n"; 

	ofstream fout(fairOutFileName);

	// space separated  values...
	for(int i = 0; i < r.size(); i++)
	{
		fout << i << " " << r[i] << "\n";
	}

	fout.close();

	return 0;
}


int DataIO :: computeFairRadius(int maxn, string fileName)
{
	for (int i=0 ; i<n ; ++i)//compute fair radiuses
	{
		r[i] = compute_r(i, n);
		// cout << i << " " << r[i] << endl;
	}

	writeFairRadiusToFile(fileName);
	cout << "Fair Radius written to file .. \n";

	return 0;
}


int DataIO :: computeWeightedFairRadius(int maxn)
{
	long double totalWeight = 0.0;

	for(int i = 0; i < corWeights.size(); i++)
	{
		totalWeight += corWeights[i];
	}
	// cout << "Total Weight = " << totalWeight << endl;
	// cout << "Total Points = " << corWeights.size() << endl;

	for (int i=0 ; i<n ; ++i)//compute fair radiuses
	{
		r[i] = compute_r_weighted(i, maxn, totalWeight);
		// cout << i << " " << r[i] << endl;
	}

	return 0;
}


int DataIO :: compute_balls(long double alpha)
{
	//computes a set of balls with respect to the fairness parameter alpha such that each point x, has a center within distance alpha*r[x]
	for(int i = 0; i < n; i++)
		covered[i] = 0;

	int center_number = 0;
	while (true)
	{
		int new_c = -1;
		long double dist = 1e100;
		for (int i=0 ; i<n ; ++i)
		{
			if (!covered[i] && r[i]<dist)
			{
				dist = r[i];
				new_c = i;
			}
		}

		if (new_c ==-1)
		{
			break;
		}
		
		center_number++;

		
		for (int i=0 ; i<n ; ++i)
		{
			// if (!covered[i] && computeEuclideanDist(point[new_c],point[i]) <= alpha*r[i]+eps)
			if (!covered[i] && compute_dist(new_c,i) <= alpha*r[i]+eps)
			{
				covered[i] = true;
 				center[i] = new_c;
			}
		}
	}
	return center_number;
}



//this function performs a binary search to find the right value of alpha for which the number of balls becomes exactly k
//this function is used in Jung et al 's algorithm for finding k centers
long double DataIO :: find_alpha()
{
	if (compute_balls(2.0)==k)
	{
		return 2.0;
	}
	
	long double b = 1.0, e=2.0;

	while (b < e-eps)
	{
		long double med = (b+e)/2.0;
		int num = compute_balls(med);
		if (num < k)
		{
			e = med;
		}
		else
		{
			b = med;
		}
	}
	if (compute_balls(b)!=k)
	{
		cerr << "------------error: not k balls in Jung's algorithm " << endl;
		b = e;
	}

	return b;
}

//this function computes the k-median cost using the array center.
long double DataIO :: compute_k_median_cost()
{
	long double ans = 0.0;
	for (int i=0; i < n; ++i)
	{
		// multiply distance by weight... -- change
		// long double d = computeEuclideanDist(point[i], point[center[i]]);
		long double d = compute_dist(i, center[i]);
		ans += d;
	}
	return ans;
}


//this function computes the weighted k-median cost using the array center.
long double DataIO :: compute_weighted_k_median_cost()
{
	long double ans = 0.0;
	for (int i=0; i < n; ++i)
	{
		// multiply distance by weight... -- change
		// long double d = computeEuclideanDist(point[i], point[center[i]]);
		long double d = compute_dist(i, center[i]);
		ans += (corWeights[i] * d);
	}
	return ans;
}


int DataIO :: computeJungEtAlAlphaAndCost()
{
	alpha = find_alpha();
	output_jung_fair = alpha;
	compute_balls(alpha);//compute critical balls
	output_jung_cost = compute_k_median_cost();	//compute cost of 

	cout << "Jung et al. alpha = " << alpha << " ----- Cost = " << output_jung_cost << endl;
	return 0;
}


int DataIO :: computeJungEtAlAlphaAndCostWeighted()
{
	alpha = find_alpha();
	output_jung_fair = alpha;
	compute_balls(alpha);//compute critical balls
	output_jung_cost = compute_weighted_k_median_cost();	//compute cost of 

	cout << "Jung et al. alpha = " << alpha << " ----- Cost = " << output_jung_cost << endl;
	return 0;
}



int DataIO :: add_center()
{//this function is used in the initialization step and adds one more center that is furthest away from current centers
	long double dist = -1.0;
	int candidate = -1;
	for (int i=0 ; i<n ; ++i)
	{
		// if (center[i]!=i && computeEuclideanDist(point[i],point[center[i]])>dist)
		if (center[i]!=i && compute_dist(i,center[i])>dist)
		{
			// dist = computeEuclideanDist(point[i],point[center[i]]);
			dist = compute_dist(i, center[i]);
			candidate = i;
		}
	}

	centers.push_back(candidate);
	iscenter[candidate] = true;

	if (whichBall[candidate]!=-1)
	{
		cnt[whichBall[candidate]]++;
	}
	for (int i=0 ; i<n ; ++i)
	{
		// if (computeEuclideanDist(point[i],point[candidate]) < computeEuclideanDist(point[i],point[center[i]]))
		if (compute_dist(i,candidate) < compute_dist(i,center[i]))
		{
			center[i] = candidate;
		}
	}

	return 0;
}


long double DataIO :: compute_cost()
{//computes k-median cost using the vectors centers
	long double ans = 0.0;
	for (int i=0 ; i<n ; ++i)
	{
		long double dis = 1e100;
		// centers is a vector of size number of clusters/centers...
		for (int j=0 ; j<centers.size() ; ++j)
		{
			// if (computeEuclideanDist(point[i], point[centers[j]]) < dis)
			if (compute_dist(i, centers[j]) < dis)
			{
				// dis = computeEuclideanDist(point[i], point[centers[j]]);
				dis = compute_dist(i, centers[j]);
			}
		}
		// take weighted cost... -- change
		ans += dis;

	}
	return ans;
}


long double DataIO :: compute_weighted_cost()
{//computes weighted k-median cost using the vectors centers
	long double ans = 0.0;
	for (int i=0 ; i<n ; ++i)
	{
		long double dis = 1e100;
		for (int j=0 ; j<centers.size() ; ++j)
		{
			// if (computeEuclideanDist(point[i], point[centers[j]]) < dis)
			if (compute_dist(i, centers[j]) < dis)
			{
				// dis = computeEuclideanDist(point[i], point[centers[j]]);
				dis = compute_dist(i, centers[j]);
			}
		}
		// take weighted cost... -- change
		// ans += dis;
		ans += (corWeights[i] * dis);
	}
	return ans;
}


//Computing the distance between the ith and jth points in the dataset
long double DataIO :: compute_dist(int i , int j)
{
	long double ans = 0.0;
	for (int dd=0 ; dd<d ; ++dd)
	{
		ans += (point[i][dd]-point[j][dd])*(point[i][dd]-point[j][dd]);
	}
	return sqrt(ans);
}


long double DataIO :: compute_dist_vectors(std::vector<long double>& point_i, std::vector<long double>& point_j)
{
	long double dist = 0.0;
	for (int dd=0 ; dd < point_i.size() ; ++dd)
	{
		dist += (point_i[dd]-point_j[dd])*(point_i[dd]-point_j[dd]);
	}
	return sqrt(dist);
}


long double DataIO :: evaluate_fairness()
{
	//computes how fair a set of centers are
	long double maxF = 0.0;
	for (int i=0 ; i<n ; ++i)
	{
		long double dis = 1e100;
		for (int j=0 ; j<centers.size() ; ++j)
		{
			// if (dis > computeEuclideanDist(point[i], point[centers[j]]))
			if (dis > compute_dist(i, centers[j]))
			{
				// dis = computeEuclideanDist(point[i], point[centers[j]]);
				dis = compute_dist(i, centers[j]);

				// if(r[i] == 0)
				// {
				// 	cout << i << " " << r[i] << " " << dis << endl;
				// }
			}
		}


		if(dis == 0 && r[i] == 0)
		{
			continue;
		}
		if (dis/r[i] > maxF)
		{
			maxF = dis/r[i];
		}
	}

	cout << "max F = " << maxF << endl;
	return maxF;
}


int DataIO :: Local_Search_initialization(long double alpha)
{//initialization for the local search algorithm
	// memset(iscenter,0,sizeof(iscenter));
	// memset(whichBall,-1,sizeof(whichBall));

	int cnum = compute_balls(3*alpha);//computes critical balls

	// cout << "cnum = " << cnum << endl;

	output_crit_balls = cnum;

	for (int i=0 ; i < n ; ++i)//setting centers to balls and other initializations
	{
		if (center[i]==i)
		{
			balls.push_back(i);
			centers.push_back(i);
			cnt.push_back(1);
			iscenter[i] = true;
		}
	}
	for (int i=0 ; i<n ; ++i)//setting value for whichBall
	{
		for (int j=0 ; j<balls.size() ; ++j)
		{
			// if (computeEuclideanDist(point[i],point[balls[j]]) < alpha*r[balls[j]])
			if (compute_dist(i,balls[j]) < alpha*r[balls[j]])
			{
				whichBall[i]=j;
			}
		}
	}
	for (int i=0 ; i<n ; ++i)
	{
		for (int j=0 ; j< balls.size() ; ++j)
		{
			// if (computeEuclideanDist(point[i], point[balls[j]]) < computeEuclideanDist(point[i], point[center[i]]))
			if (compute_dist(i, balls[j]) < compute_dist(i, center[i]))
				center[i] = balls[j];
		}
	}

	for (int i=cnum ; i<k ; ++i)
	{
		add_center();
	}

	return 0;
}


long double DataIO :: Local_Search(long double alpha)
{
	//local search algorithm
	Local_Search_initialization(alpha);
	output_init_cost = compute_cost();
	output_init_fair = evaluate_fairness();

	cout << "init ls cost = " << output_init_cost << "init ls fair = " << output_init_fair << endl;

	bool improved = true;
	while (improved)
	{
		improved = false;
		for (int i=0 ; i<centers.size() ; ++i)
		{
			for (int j=0 ; j<n ; ++j)
			{
				if (!iscenter[j])
				{
					int b = whichBall[centers[i]];
					if (b==-1 || cnt[b]>1 || whichBall[j]==b)
					{//checks for feasibility
						long double prev_cost = compute_cost();
						int prev_c = centers[i];
						centers[i] = j;
						long double cc = compute_cost();
						if (cc<prev_cost*0.99)
						{//checks if the cost is improved
							iscenter[j] = true;
							iscenter[prev_c] = false;
							
							if (whichBall[prev_c]!=-1)
							{
								cnt[whichBall[prev_c]]--;
							}
							
							if (whichBall[j]!=-1)
							{
								cnt[whichBall[j]]++;
							}
							improved = true;
							std::cout << "Improved cost -- " << cc << std::endl;
						}
						else
						{
							centers[i] = prev_c;
						}
					}
				}
			}
		}
	}
	return compute_cost();
}


long double DataIO :: Local_Search_Weighted(long double alpha)
{
	//local search algorithm
	Local_Search_initialization(alpha);
	output_init_cost = compute_weighted_cost();
	output_init_fair = evaluate_fairness();

	cout << "init ls cost = " << output_init_cost << "init ls fair = " << output_init_fair << endl;

	bool improved = true;
	while (improved)
	{
		improved = false;
		for (int i=0 ; i<centers.size() ; ++i)
		{
			for (int j=0 ; j<n ; ++j)
			{
				if (!iscenter[j])
				{
					int b = whichBall[centers[i]];
					if (b==-1 || cnt[b]>1 || whichBall[j]==b)
					{//checks for feasibility
						long double prev_cost = compute_weighted_cost();
						int prev_c = centers[i];
						centers[i] = j;
						long double cc = compute_weighted_cost();
						if (cc<prev_cost*0.99)
						{//checks if the cost is improved
							iscenter[j] = true;
							iscenter[prev_c] = false;
							
							if (whichBall[prev_c]!=-1)
							{
								cnt[whichBall[prev_c]]--;
							}
							
							if (whichBall[j]!=-1)
							{
								cnt[whichBall[j]]++;
							}
							improved = true;
							std::cout << "Improved cost -- " << cc << std::endl;
						}
						else
						{
							centers[i] = prev_c;
						}
					}
				}
			}
		}
	}
	return compute_weighted_cost();
}


int DataIO :: storeCenters()
{
	cout << "Printing centers -- " << centers.size() << endl;
	for(int i = 0; i < centers.size(); i++)
	{
		cout << i << " " << centers[i] << " " <<  printVec(point[centers[i]]) << endl;
	}

	unordered_map<int, bool> uniqueCenterIds; 

	// printVec(point[center[0]]);
	for(int i = 0; i < center.size(); i++)
	{
		uniqueCenterIds[center[i]] = 1;
	}
	cout << "Got the unique centers...\n"; 
	unordered_map<int, bool> :: iterator it;

	for(it = uniqueCenterIds.begin(); it != uniqueCenterIds.end(); ++it)
	{

	// for(int i = 0; i < center.size(); i++)
	// {
		centerVectors.push_back(it->first);
	}

	cout << "Printing center vectors -- " << centerVectors.size() << endl;
	for(int i = 0; i < centerVectors.size(); i++)
	{
		cout << i << " " << centerVectors[i] << " "<< printVec(point[centerVectors[i]]) << endl;
	}
	return 0;
}

int DataIO :: computeMVIndFair()
{
	cout << "MV Ind Fairness --- \n";
	output_ls_cost = Local_Search(alpha);
	cout << "Evaluating Fairness...\n";
	output_ls_fair = evaluate_fairness();

	// storeCenters();
	// cout << "stored centers\n";
	return 0;
}

int DataIO :: computeMVIndFairWeighted()
{
	cout << "MV Ind Fairness --- \n";
	output_ls_cost = Local_Search_Weighted(alpha);
	cout << "Evaluating Fairness...\n";
	output_ls_fair = evaluate_fairness();

	// storeCenters();
	// cout << "stored centers\n";
	return 0;
}


long double DataIO :: computeCostForPoint(vector<long double> eachPoint)
{
	// long double minCost = computeEuclideanDist(eachPoint, point[centerVectors[0]]);
	long double minCost = compute_dist_vectors(eachPoint, point[centers[0]]);

	for(int i = 1; i < centers.size(); i++)
	{
		// long double thisCenterCost = computeEuclideanDist(eachPoint, point[centerVectors[i]]);
		long double thisCenterCost = compute_dist_vectors(eachPoint, point[centers[i]]);
		if (minCost > thisCenterCost)
		{
			minCost = thisCenterCost;
		} 
	}

	return minCost;
}


long double DataIO :: readFullFileAndComputeCost(string fileName)
{

	long double totalCost = 0.0;
	stringstream ss;
	string line;

	ifstream fin(fileName);
	// space separated  values...

	while(getline(fin, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		vector<long double> tempVec;
		long double tempval;
		while(ss >> tempval)
		{
			tempVec.push_back(tempval);
		}

		long double costOnPoint = computeCostForPoint(tempVec);		// this is just the euclidean distance...

		totalCost += costOnPoint;

		distOfEachPointFullData.push_back(costOnPoint);
	}

	// cout << "Size - distOfEachPointFullData = " << distOfEachPointFullData.size() << endl;

	fin.close();
	// printVecVec(point);
	return totalCost;
}


int DataIO :: computeCostOnFullDataUsingCoresetCenters(string fileName)
{
	costOnFullData = readFullFileAndComputeCost(fileName);
	 // = std::accumulate(localDistOfEachPoint.begin(), localDistOfEachPoint.end(), 0.0);
	cout << "Got the cost on full data --- " << costOnFullData << endl;
	return 0;
}

int DataIO :: computeFairnessOnFullDataUsingCoresetCenters(string fairRadiusFileName)
{

	long double maxFairness = 0, fairnessForPoint;
	stringstream ss;
	string line;

	ifstream fin(fairRadiusFileName);
	// space separated  values...
	int i = 0;
	while(getline(fin, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		int pointId;
		long double pointFairRadius;
		
		ss >> pointId >> pointFairRadius;
		
		fairnessForPoint = distOfEachPointFullData[i]/pointFairRadius;

		if(maxFairness < fairnessForPoint)
		{
			maxFairness = fairnessForPoint;
		}

		i += 1;
	}

	fin.close();

	fairnessOnFullData = maxFairness;
	if(maxFairness == 0)
	{
		cout << "Here it is...\n";
		exit(0);
	}
	cout << "Fairness evaluated --" << maxFairness << endl;
	return 0;
}


int DataIO :: generate_output(string fileName, double computeTime, int flag)
{
	vector<string> pathFileName = getPathAndFileName(fileName);

	string strk = std::to_string(k);

	string outFileName = "output/out_" + strk + "_" + pathFileName[1];
	outFileName = pathFileName[0] + "/" + outFileName;

	// fileName = "./allOutput/out_" + fileName;

	cout << "writing output to -- " << outFileName << endl;

	ofstream fout(outFileName);

	if(flag == 0)
	{
		fout << "n,d,k,output_crit_balls,output_jung_fair,output_jung_cost,output_init_fair,output_init_cost,output_ls_fair,output_ls_cost,computeTime\n";
	}
	else
	{
		fout << "n,d,k,output_crit_balls,output_jung_fair,output_jung_cost,output_init_fair,output_init_cost,output_ls_fair,output_ls_cost,fairnessOnFullData,costOnFullData,computeTime\n";
	}

	fout << n << ",";
	fout << d << ",";
	fout << k << ",";
	
	fout << output_crit_balls << ",";
	
	fout << output_jung_fair << ",";
	fout << setprecision(20) << output_jung_cost << ",";
	
	fout << output_init_fair << ",";
	fout << setprecision(20) << output_init_cost << ",";

	fout << output_ls_fair << ",";
	fout << setprecision(20) << output_ls_cost << ",";

	if(flag)
	{
		fout << fairnessOnFullData << ",";
		fout << setprecision(20) << costOnFullData << ",";
	}

	fout << setprecision(20) << computeTime << "\n";

	fout.close();

	return 0;
}