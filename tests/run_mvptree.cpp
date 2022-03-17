/**
 * g++ -omvptree run_mvptree.cpp -g -O2 -Wall -Wno-unused-variable -I include/
 **/
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>
#include <chrono>
#include <ratio>
#include <mvptree/mvptree.hpp>
#include <mvptree/datapoint.hpp>

using namespace std;


static long long m_id = 1;
static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);

static const int n_centers = 10;

const int bf = 2;   //branchfactor
const int pl = 2;   // pathlength
const int lc = 2500; // leafcap
const int lpn = 2;  // levelspernode
const int fo = 4; //fanout bf^lpn
const int ns = 2; //numsplits bf^(lpn-1)

const int n_runs = 10;

const int n = 100000;
const int n_iters = 5;

const int cluster_size = 10;
const int radius = 10;


struct PerfMetrics {
	double build_ops;
	double build_time;
	double query_ops;
	double query_time;
	PerfMetrics(){};
	PerfMetrics(const PerfMetrics &other){
		build_ops = other.build_ops;
		build_time = other.build_time;
		query_ops = other.query_ops;
		query_time = other.query_time;
	}
};

int generate_data(vector<DataPoint<pl>*> &points, int N){
	for (int i=0;i<N;i++){
		uint64_t value = m_distrib(m_gen);
		DataPoint<pl> *pnt = new H64DataPoint<pl>(m_id++, value);
		points.push_back(pnt);
	}
	return 0;
}

int generate_cluster(vector<DataPoint<pl>*> &points, uint64_t center, int N, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;

	DataPoint<pl> *mid = new H64DataPoint<pl>(m_id++, center);
	points.push_back(mid);
	for (int i=0;i < N-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		DataPoint<pl> *pnt = new H64DataPoint<pl>(m_id++, val);
		points.push_back(pnt);
	}
	return N;
}

void do_run(MVPTree<bf,pl,lc,lpn,fo,ns> &mvptree, vector<PerfMetrics> &metrics){

	DataPoint<pl>::ResetCount();

	chrono::duration<double> total(0);
	for (int i=0;i < n_iters;i++){
		vector<DataPoint<pl>*> points;
		generate_data(points, n);
		
		auto s = chrono::steady_clock::now();
		mvptree.Add(points);
		auto e = chrono::steady_clock::now();
		total += (e - s);
	}

	mvptree.Sync();
	int sz = mvptree.Size();

	PerfMetrics avgs;
	avgs.build_ops  = ((double)DataPoint<pl>::n_ops/(double)sz)*100.0;
	avgs.build_time = total.count();
	
	cout << "build: "  << setw(10) << setprecision(6) << avgs.build_ops << "% opers - ("
		 << setw(10) << setprecision(6) << avgs.build_time << " seconds)";

	uint64_t centers[n_centers];

	for (int i=0;i < n_centers;i++){
		centers[i] = m_distrib(m_gen);
		vector<DataPoint<pl>*> points;
		generate_cluster(points, centers[i], cluster_size, radius);
		mvptree.Add(points);
	}

	mvptree.Sync();
	
	sz = mvptree.Size();


	DataPoint<pl>::ResetCount();
	chrono::duration<double, milli> total_query(0);
	for (int i=0;i < n_centers;i++){
		H64DataPoint<pl> target(0, centers[i]);
		auto s = chrono::steady_clock::now();
		vector<DataPoint<pl>*> results = mvptree.Query(target, radius);
 		auto e = chrono::steady_clock::now();
		total_query += (e - s);
	}

	avgs.query_ops = ((double)DataPoint<pl>::n_ops/(double)sz/(double)n_centers)*100.0;
	avgs.query_time = (double)total_query.count()/(double)n_centers;
	
	cout << setw(10) << right << "query: "
		 << setw(10)  << setprecision(6) << avgs.query_ops << "% distance operations - ("
		 << setw(10) << setprecision(6)  << avgs.query_time << " millisecs)" <<  endl;

	metrics.push_back(avgs);
	mvptree.Clear();
	return;
}


int main(int argc, char **argv){


	
	MVPTree<bf,pl,lc,lpn,fo,ns> mvptree;

	cout << "MVPTree<bf,pl,lc,lpn,fo,ns>  ";
	cout << "<" << bf << "," << pl << "," << lc << "," << lpn << "," << fo << "," << ns << ">" << endl << endl;

	vector<PerfMetrics> metrics;
	for (int j=0;j < n_runs;j++){
		cout <<  "(" << j+1 << ")";
		do_run(mvptree, metrics);
	}
	cout << endl << endl;


	double avg_build_ops = 0;
	double avg_build_time = 0;
	double avg_query_ops = 0;
	double avg_query_time = 0;

	for (PerfMetrics curr : metrics){
		avg_build_ops += curr.build_ops/n_runs;
		avg_build_time += curr.build_time/n_runs;
		avg_query_ops += curr.query_ops/n_runs;
		avg_query_time += curr.query_time/n_runs;
	}
	
	cout << endl << "avg: build ops: " << avg_build_ops << "% ops"
		 << " build time: " << avg_build_time << " seconds" << endl
		 << " avg query ops: " << avg_query_ops << "% ops" 
		 << " avg query time: " << avg_query_time << " millisecs" << endl;

	
	cout << "Done." << endl;
	return 0;
}
