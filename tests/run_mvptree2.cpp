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

const int n_runs = 4;

const int npoints = 50000;
const int n_iters = 10;
const int nclusters = 10;
const int cluster_size = 10;
const double radius = 0.04;
const double diff = sqrt(pow(radius, 2.0)/16.0);

const int BF = 2;   //branchfactor
const int PL = 8;   // pathlength
const int LC = 250; // leafcap
const int LPN = 10;  // levelspernode
const int FO = 1024; //fanout bf^lpn
const int NS = 512; //numsplits bf^(lpn-1)

static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_real_distribution<double> m_distrib(-1.0, 1.0);
static 	uniform_real_distribution<double> radius_distr(-diff, diff);


struct perfmetric {
	double avg_build_ops;
	double avg_build_time;
	double avg_query_ops;
	double avg_query_time;
};

int generate_point(double v[]){

	for (int i=0;i < 16;i++) v[i] = m_distrib(m_gen);

	return 0;
}

int generate_data(vector<datapoint_t<VectorKeyObject,PL>> &points, int n){
	double v[16];
	for (int i=0;i < n;i++){
		generate_point(v);
		points.push_back({ m_id++, VectorKeyObject(v) });
	}
	return 0;
}

int generate_cluster(vector<datapoint_t<VectorKeyObject,PL>> &points, double center[], int n){

	points.push_back({ g_id++, VectorKeyObject(center) });

	double v[16];
	for (int i=0;i < n-1;i++){
		memcpy(v, center, 16*sizeof(double));
		for (int j=0;j < 16;j++) v[j] += radius_distr(m_gen);
		points.push_back({ g_id++, VectorKeyObject(v) });
	}
	return 0;
}

void do_run(int index, vector<struct perfmetric> &metrics){
	m_id = 1;
	g_id = 100000000;

	MVPTree<VectorKeyObject,BF,PL,LC,LPN,FO,NS> tree;

	int sz;
	chrono::duration<double> total(0);
	datapoint_t<VectorKeyObject,PL>::n_build_ops = 0;
	for (int i=0;i < n_iters;i++){
		vector<datapoint_t<VectorKeyObject,PL>> points;
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
	
	m.avg_build_ops = 100.0*((double)datapoint_t<VectorKeyObject,PL>::n_build_ops/(double)sz);
	m.avg_build_time = total.count();
	
	cout << "(" << index << ") build tree: " << setw(10) << setprecision(6) << m.avg_build_ops << "% opers "
		 << setw(10) << setprecision(6) << m.avg_build_time << " secs ";

	double centers[nclusters][16];
	for (int i=0;i < nclusters;i++){
		generate_point(centers[i]);

		vector<datapoint_t<VectorKeyObject,PL>> cluster;
		generate_cluster(cluster, centers[i], cluster_size);
		assert(cluster.size() == cluster_size);
		
		tree.Add(cluster);
		sz = tree.Size();

		assert(sz == npoints*n_iters + cluster_size*(i+1));
	}

	tree.Sync();

	datapoint_t<VectorKeyObject,PL>::n_query_ops = 0;
	chrono::duration<double, milli> querytime(0);
	for (int i=0;i < nclusters;i++){
		VectorKeyObject target(centers[i]);

		auto s = chrono::steady_clock::now();
		vector<item_t<VectorKeyObject>> results = tree.Query(target, radius);
		auto e = chrono::steady_clock::now();
		querytime += (e - s);

		int nresults = (int)results.size();
		assert(nresults >= cluster_size);
	}

	m.avg_query_ops = 100.0*((double)datapoint_t<VectorKeyObject,PL>::n_query_ops/(double)sz/(double)nclusters);
	m.avg_query_time = (double)querytime.count()/(double)nclusters;

	cout << " query ops " << dec << setprecision(6) << m.avg_query_ops << "% opers   " 
		 << "query time: " << dec << setprecision(6) << m.avg_query_time << " millisecs" << endl;

	metrics.push_back(m);

	tree.Clear();
	assert(tree.Size() == 0);
}


int main(int argc, char **argv){

	MVPTree<VectorKeyObject,BF,PL,LC,LPN,FO,NS> tree;

	cout << "MVPTree Tests" << endl << endl;;

	cout << "tree parameters: " << endl;
	cout << "--------------- " << endl;
	cout << "branchfactor, bf = " << BF << endl;
	cout << "path length, pl = " << PL << endl;
	cout << "leaf capacity, lc = " << LC << endl;
	cout << "levels per node, lpn = " << LPN << endl;
	cout << "internal node fanout, fo = " << FO << endl << endl;
	
	cout << "10-dim real valued object vectors in euclidean metric space" << endl;
	cout << "dataset size: " << npoints*n_iters << " datapoints" << endl;
	cout << "no. test runs averaged over: " << n_runs << endl;
	
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
