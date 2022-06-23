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
#include <mvptree/key.hpp>

using namespace std;
using namespace mvp;

static long long m_id = 1;
static long long g_id = 100000000;

static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);


const int BF = 2;   //branchfactor
const int PL = 8;   // pathlength
const int LC = 100; // leafcap
const int LPN = 4;  // levelspernode
const int FO = 16; //fanout bf^lpn
const int NS = 8; //numsplits bf^(lpn-1)


struct perfmetric {
	double avg_build_ops;
	double avg_build_time;
	double avg_query_ops;
	double avg_query_time;
	double avg_memory;
};

int generate_data(vector<datapoint_t<H64KeyObject,PL>> &points, int n){
	for (int i=0;i < n;i++){
		uint64_t value = m_distrib(m_gen);
		points.push_back({ m_id++, H64KeyObject(value) });
	}
	return 0;
}

int generate_cluster(vector<datapoint_t<H64KeyObject,PL>> &points, uint64_t center, int n, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;
	points.push_back({ g_id++, H64KeyObject(center) });
	for (int i=0;i < n-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		points.push_back({ g_id++, H64KeyObject(val) });
	}
	return 0;
}

void do_run(const int index, const int n_points, const int n_clusters,
			const int clustersize, const int radius, vector<struct perfmetric> &metrics){
	m_id = 1;
	g_id = 100000000;

	MVPTree<H64KeyObject,BF,PL,LC,LPN,FO,NS> tree;

	int sz = 0;
	const int lotsz = 1000;
	chrono::duration<double, micro> total(0);
	datapoint_t<H64KeyObject,PL>::n_build_ops = 0;

	while (sz < n_points){
		vector<datapoint_t<H64KeyObject,PL>> points;
		generate_data(points, lotsz);
		auto s = chrono::steady_clock::now();
		tree.Add(points);
		auto e = chrono::steady_clock::now();
		total += (e - s);
		sz = tree.Size();
	}
	tree.Sync();
	assert(sz == n_points);

	struct perfmetric m;
	m.avg_build_ops = 100.0*((double)datapoint_t<H64KeyObject,PL>::n_build_ops/(double)sz);
	m.avg_build_time = total.count()/(double)sz;
	m.avg_memory = tree.MemoryUsage();
	
	cout << "(" << index << ") build tree: " << setw(10) << setprecision(6) << m.avg_build_ops << "% opers "
		 << setw(10) << setprecision(6) << m.avg_build_time << " microsecs ";

	uint64_t centers[n_clusters];
	for (int i=0;i < n_clusters;i++){
		centers[i] = m_distrib(m_gen);

		vector<datapoint_t<H64KeyObject,PL>> cluster;
		generate_cluster(cluster, centers[i], clustersize, radius);
		assert(cluster.size() == clustersize);
		
		tree.Add(cluster);
		sz = tree.Size();
		assert(sz == npoints + clustersize*(i+1));
	}

	tree.Sync();

	datapoint_t<H64KeyObject,PL>::n_query_ops = 0;
	chrono::duration<double, milli> querytime(0);
	for (int i=0;i < n_clusters;i++){
		H64KeyObject target(centers[i]);

		auto s = chrono::steady_clock::now();
		vector<item_t<H64KeyObject>> results = tree.Query(target, radius);
		auto e = chrono::steady_clock::now();
		querytime += (e - s);

		int nresults = (int)results.size();
		assert(nresults >= clustersize);
	}

	m.avg_query_ops = 100.0*(double)datapoint_t<H64KeyObject,PL>::n_query_ops/(double)n_clusters/(double)sz;
	m.avg_query_time = (double)querytime.count()/(double)n_clusters;

	cout << " query ops " << dec << setprecision(6) << m.avg_query_ops << "% opers   " 
		 << "query time: " << dec << setprecision(6) << m.avg_query_time << " millisecs" << endl;

	metrics.push_back(m);

	tree.Clear();
	assert(tree.Size() == 0);
}

void do_experiment(const int n, const int n_runs, const int n_points,
				   const int n_clusters, const int clustersize, const int radius){

	cout << "-------------Experiment = " << n << "-----------------" << endl;
	cout << "size: " << n_points << endl;
	cout << n_clusters << " clusters of " << clustersize << " with radius = " << radius << endl;
	cout << "no. runs: " << n_runs << endl;

	vector<perfmetric> metrics;
	for (int i=0;i < n_runs;i++){
		do_run(i, n_points, n_clusters, clustersize, radius, metrics);
	}

	double avg_build_ops = 0;
	double avg_build_time = 0;
	double avg_query_ops = 0;
	double avg_query_time = 0;
	double avg_memory = 0;
	for (struct perfmetric &m : metrics){
		avg_build_ops += m.avg_build_ops/n_runs;
		avg_build_time += m.avg_build_time/n_runs;
		avg_query_ops += m.avg_query_ops/n_runs;
		avg_query_time += m.avg_query_time/n_runs;
		avg_memory += m.avg_memory/n_runs;
	}

	cout << "avg build:  " << avg_build_ops << "% opers " << avg_build_time << " microseconds" << endl;
	cout << "avg query:  " << avg_query_ops << "% opers " << avg_query_time << " milliseconds" << endl;
	cout << "avg memory: " << avg_memory/1000000.0 << " MB" << endl;
	cout << "--------------------------------------------------" << endl << endl;
	
}

int main(int argc, char **argv){

	MVPTree<H64KeyObject,BF,PL,LC,LPN,FO,NS> tree;

	cout << "MVPTree Tests" << endl << endl;;

	cout << "tree parameters: " << endl;
	cout << "--------------- " << endl;
	cout << "branchfactor, bf = " << BF << endl;
	cout << "path length, pl = " << PL << endl;
	cout << "leaf capacity, lc = " << LC << endl;
	cout << "levels per node, lpn = " << LPN << endl << endl;

	const int n_runs = 5;
	const int n_experiments = 7;
	const int npoints[n_experiments] = { 100000, 200000, 400000, 800000,
										 1000000, 2000000, 4000000};
	const int nclusters = 10;
	const int clustersize = 10;
	const int radius = 10;

	for (int i=0;i < n_experiments;i++){
		do_experiment(i+1, n_runs, npoints[i], nclusters, clustersize, radius);
	}
	
	cout << "Experiments with various radius settings" << endl;
	
	const int N = 4000000;
	const int n_rad = 7;
	const int rad[n_rad] = { 0, 2, 4, 6, 8, 10, 12 };

	for (int i=0;i < n_rad;i++){
		do_experiment(i+1, n_runs, N, nclusters, clustersize, rad[i]);
	}
	
	return 0;
}
