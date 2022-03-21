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
#include <cassert>
#include <ratio>
#include <mvptree/mvptree.hpp>
#include <mvptree/datapoint.hpp>

using namespace std;


static long long m_id = 1;
static long long g_id = 100000000;

static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);


const int n_runs = 4;

const int npoints = 100000;
const int n_iters = 10;
const int nclusters = 10;
const int cluster_size = 10;
const int radius = 10;

const int BF = 2;   //branchfactor
const int PL = 10;   // pathlength
const int LC = 3000; // leafcap
const int LPN = 10;  // levelspernode
const int FO = 1024; //fanout bf^lpn
const int NS = 512; //numsplits bf^(lpn-1)


struct perfmetric {
	double avg_build_ops;
	double avg_build_time;
	double avg_query_ops;
	double avg_query_time;
};


int generate_data(vector<DataPoint<PL>*> &points, int n){
	for (int i=0;i < n;i++){
		uint64_t value = m_distrib(m_gen);
		DataPoint<PL> *pnt = new H64DataPoint<PL>(m_id++, value);
		points.push_back(pnt);
	}
	return 0;
}

int generate_cluster(vector<DataPoint<PL>*> &points, uint64_t center, int n, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;

	DataPoint<PL> *mid = new H64DataPoint<PL>(g_id++, center);
	points.push_back(mid);
	for (int i=0;i < n-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		DataPoint<PL> *pnt = new H64DataPoint<PL>(g_id++, val);
		points.push_back(pnt);
	}
	return 0;
}

void do_run(int index, vector<struct perfmetric> &metrics){
	m_id = 1;
	g_id = 100000000;

	MVPTree<BF,PL,LC,LPN,FO,NS> tree;

	int sz;
	chrono::duration<double> total(0);
	DataPoint<PL>::n_ops = 0;
	for (int i=0;i < n_iters;i++){
		vector<DataPoint<PL>*> points;
		generate_data(points, npoints);
	
		auto s = chrono::steady_clock::now();
		tree.Add(points);
		auto e = chrono::steady_clock::now();
		total += (e - s);
		sz = tree.Size();

		assert(sz == npoints*(i+1));
	}

	tree.Sync();

	struct perfmetric m;
	
	m.avg_build_ops = 100.0*((double)DataPoint<PL>::n_ops/(double)sz);
	m.avg_build_time = total.count();
	
	cout << "(" << index << ") build tree: " << setw(10) << setprecision(6) << m.avg_build_ops << "% opers "
		 << setw(10) << setprecision(6) << m.avg_build_time << " secs ";

	uint64_t centers[nclusters];
	for (int i=0;i < nclusters;i++){
		centers[i] = m_distrib(m_gen);

		vector<DataPoint<PL>*> cluster;
		generate_cluster(cluster, centers[i], cluster_size, radius);
		assert(cluster.size() == cluster_size);
		
		tree.Add(cluster);
		sz = tree.Size();

		assert(sz == npoints*n_iters + cluster_size*(i+1));
	}

	tree.Sync();

	DataPoint<PL>::n_ops = 0;
	chrono::duration<double, milli> querytime(0);
	for (int i=0;i < nclusters;i++){
		H64DataPoint<PL> target(0, centers[i]);

		auto s = chrono::steady_clock::now();
		vector<DataPoint<PL>*> results = tree.Query(target, radius);
		auto e = chrono::steady_clock::now();
		querytime += (e - s);

		int nresults = (int)results.size();
		assert(nresults >= cluster_size);
	}

	m.avg_query_ops = 100.0*((double)DataPoint<PL>::n_ops/(double)sz/(double)nclusters);
	m.avg_query_time = (double)querytime.count()/(double)nclusters;

	cout << " query ops " << dec << setprecision(6) << m.avg_query_ops << "% opers   " 
		 << "query time: " << dec << setprecision(6) << m.avg_query_time << " millisecs" << endl;

	metrics.push_back(m);

	tree.Clear();
	assert(tree.Size() == 0);
}


int main(int argc, char **argv){

	MVPTree<BF,PL,LC,LPN,FO,NS> tree;

	cout << "    MVPTree<bf,pl,lc,lpn,fo,ns>  ";
	cout << "           <" << BF << "," << PL << "," << LC << "," << LPN << "," << FO << "," << NS << ">" << endl;

	vector<perfmetric> metrics;
	for (int i=0;i < n_runs;i++){
		do_run(i, metrics);
	}

	double avg_build_ops = 0;
	double avg_build_time = 0;
	double avg_query_ops = 0;
	double avg_query_time = 0;
	for (struct perfmetric &m : metrics){
		avg_build_ops += m.avg_build_ops/n_runs;
		avg_build_time += m.avg_build_time/n_runs;
		avg_query_ops += m.avg_query_ops/n_runs;
		avg_query_time += m.avg_query_time/n_runs;
	}

	cout << "no. runs: " << metrics.size() << endl;
	
	cout << "avg build:  " << avg_build_ops << "% opers " << avg_build_time << " seconds" << endl;
	cout << "avg query:  " << avg_query_ops << "% opers " << avg_query_time << " milliseconds" << endl;

	
	cout << endl << "Done." << endl;
	return 0;
}
