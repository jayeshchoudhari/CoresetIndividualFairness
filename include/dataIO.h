#ifndef DATAIO_H
#define DATAIO_H

class DataIO
{
	private:
	public:
		int n;		// Number of datapoints
		int d;		// Dimension of the datapoints
		int k; 		// Number of Centers

		long double eps;

		std::vector <std::vector <long double> > point;
		std::vector<long double> r;			//fair radius of each point

		long double alpha;

		std::vector<bool> covered;//shows if a point is covered by current centers
		std::vector<int> center;//shows the corresponding center for each point

		std::vector <int> balls; //keeps the critical balls used by our algorithm
		std::vector <int> centers; // the set of centers found by our algorithm
		std::vector <int> cnt; //keeps track of the number of centers in each ball


		long double output_init_fair, output_init_cost; //keeps the cost and fairness of the solution right after the initialization step
		long double output_jung_fair, output_jung_cost; //fairness and cost of Jung et al.'s algorithm
		long double output_ls_fair, output_ls_cost;//fairness and cost of local search

		std::vector<int> whichBall;//keeps for each point in which critical ball they exist or -1 otherwise
		std::vector<bool> iscenter;//keeps for each point if it is in the current set of centers
		int output_crit_balls;//for outputting the number of critical balls

		DataIO(int n, int d, int k, std::string inputFileName, int maxn, int maxd, long double eps);

		int readInput(std::string filename);

		int computeFairRadius(int maxn);
		long double compute_r(int i, int maxn);

		int computeJungEtAlAlphaAndCost();
		long double find_alpha();
		int compute_balls(long double alpha);

		long double compute_k_median_cost();

		int computeMVIndFair();
		int add_center();
		long double evaluate_fairness();
		long double compute_cost();
		int Local_Search_initialization(long double alpha);
		long double Local_Search(long double alpha);

		int generate_output();

};


#endif //DATAIO_H