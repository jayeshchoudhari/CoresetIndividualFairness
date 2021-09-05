#include "namespace.h"
#include "dataIO.h"
#include "util.h"

using namespace std;

DataIO :: DataIO(int n, int d, int k, string inputFileName, int maxn, int maxd, long double eps)
{
	n = n;
	d = d;
	k = k;
	eps = eps;

	point.resize(maxn, vector<long double>(maxd));

	r.resize(maxn);

	covered.resize(maxn, 0);	//shows if a point is covered by current centers
	center.resize(maxn, -1);	//shows the corresponding center for each point -- default set to -1

	whichBall.resize(maxn, -1);	//keeps for each point in which critical ball they exist or -1 otherwise
	iscenter.resize(maxn, 0);	//keeps for each point if it is in the current set of centers

	readInput(inputFileName);
}

//Reading the input points. Each point comes in a separate line and the coordinates are separated by comma.
int DataIO :: readInput(string filename)
{
	// ifstream fin("input.csv");
	ifstream fin(filename);
	for (int i=0 ; i<n ; ++i)
	{
		for (int j=0 ; j<d ; ++j)
		{
			fin >> point[i][j];
			char c;
			if (j<d-1)
				fin >> c;		// for the comma...
		}
	}
	return 0;
}

long double DataIO :: compute_r(int i, int maxn)
{
	//computes the fair radius for each point
	long double temp[maxn];

	for (int j=0 ; j<n ; ++j)
	{
		temp[j] = computeEuclideanDist(point[i], point[j]);
	}

	sort(temp,temp+n);

	// have to consider weights... 
	// take sum of first l points that weigh to n-1/k -- change
	return temp[(n-1)/k];
}


int DataIO :: computeFairRadius(int maxn)
{
	for (int i=0 ; i<n ; ++i)//compute fair radiuses
		r[i] = compute_r(i, maxn);

	return 0;
}


int DataIO :: compute_balls(long double alpha)
{
	//computes a set of balls with respect to the fairness parameter alpha such that each point x, has a center within distance alpha*r[x]
	// memset(covered, 0, sizeof(covered));

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
			if (!covered[i] && computeEuclideanDist(point[new_c],point[i]) <= alpha*r[i]+eps)
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
		long double d = computeEuclideanDist(point[i], point[center[i]]);
		ans += d;
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



int DataIO :: add_center()
{//this function is used in the initialization step and adds one more center that is furthest away from current centers
	long double dist = -1.0;
	int candidate = -1;
	for (int i=0 ; i<n ; ++i)
	{
		if (center[i]!=i && computeEuclideanDist(point[i],point[center[i]])>dist)
		{
			dist = computeEuclideanDist(point[i],point[center[i]]);
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
		if (computeEuclideanDist(point[i],point[candidate])<computeEuclideanDist(point[i],point[center[i]]))
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
		for (int j=0 ; j<centers.size() ; ++j)
		{
			if (computeEuclideanDist(point[i], point[centers[j]]) < dis)
			{
				dis = computeEuclideanDist(point[i], point[centers[j]]);
			}
		}
		// take weighted cost... -- change
		ans += dis;
	}
	return ans;
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
			if (dis > computeEuclideanDist(point[i], point[centers[j]]))
			{
				dis = computeEuclideanDist(point[i], point[centers[j]]);
			}
		}
		if (dis/r[i] > maxF)
		{
			maxF = dis/r[i];
		}
	}
	return maxF;
}


int DataIO :: Local_Search_initialization(long double alpha)
{//initialization for the local search algorithm
	// memset(iscenter,0,sizeof(iscenter));
	// memset(whichBall,-1,sizeof(whichBall));

	int cnum = compute_balls(3*alpha);//computes critical balls

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
			if (computeEuclideanDist(point[i],point[balls[j]]) < alpha*r[balls[j]])
			{
				whichBall[i]=j;
			}
		}
	}
	for (int i=0 ; i<n ; ++i)
	{
		for (int j=0 ; j< balls.size() ; ++j)
		{
			if (computeEuclideanDist(point[i], point[balls[j]]) < computeEuclideanDist(point[i], point[center[i]]))
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
						if (compute_cost()<prev_cost*0.99)
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


int DataIO :: computeMVIndFair()
{
	output_ls_cost = Local_Search(alpha);
	output_ls_fair = evaluate_fairness();
}

int DataIO :: generate_output()
{
	ofstream fout("output.csv", ofstream::app);
	fout << n << ",";
	fout << d << ",";
	fout << k << ",";
	
	fout << output_crit_balls << ",";
	
	fout << output_jung_fair << ",";
	fout << output_jung_cost << ",";
	fout << output_init_fair << ",";
	fout << output_init_cost << ",";
	fout << output_ls_fair << ",";
	fout << output_ls_cost << endl;
	fout.close();

	return 0;
}