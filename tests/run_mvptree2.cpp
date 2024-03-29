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
#include <cstring>
#include <mvptree/mvptree.hpp>
#include <mvptree/datapoint.hpp>
#include <mvptree/key.hpp>

using namespace std;
using namespace mvp;

static long long m_id = 1;
static long long g_id = 100000000;


const int BF = 2;   //branchfactor
const int PL = 8;   // pathlength
const int LC = 50; // leafcap
const int LPN = 4;  // levelspernode
const int FO = 16; //fanout bf^lpn
const int NS = 8; //numsplits bf^(lpn-1)

static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_real_distribution<double> m_distrib(-1.0, 1.0);


struct perfmetric {
	double avg_build_ops;
	double avg_build_time;
	double avg_query_ops;
	double avg_query_time;
	double avg_memory;
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

int generate_cluster(vector<datapoint_t<VectorKeyObject,PL>> &points, double center[], int n, double radius){

	double diff = sqrt(pow(radius, 2.0)/16.0);
 	uniform_real_distribution<double> radius_distr(-diff, diff);

	points.push_back({ g_id++, VectorKeyObject(center) });

	double v[16];
	for (int i=0;i < n-1;i++){
		memcpy(v, center, 16*sizeof(double));
		if (radius > 0){
			for (int j=0;j < 16;j++) v[j] += radius_distr(m_gen);
		}
		points.push_back({ g_id++, VectorKeyObject(v) });
	}
	return 0;
}

void do_run(const int index, const int n_points, const int n_clusters,
			const int clustersize,  const double radius, vector<struct perfmetric> &metrics){

	m_id = 1;
	g_id = 100000000;

	MVPTree<VectorKeyObject,BF,PL,LC,LPN,FO,NS> tree;

	size_t sz = 0;
	chrono::duration<double,micro> total(0);
	datapoint_t<VectorKeyObject,PL>::n_build_ops = 0;

	const int lot_size = 100;
	int i = 1;
	while ((int)sz < n_points){
		vector<datapoint_t<VectorKeyObject,PL>> points;
		generate_data(points, lot_size);
	
		auto s = chrono::steady_clock::now();
		tree.Add(points);
		auto e = chrono::steady_clock::now();
		total += (e - s);
		sz = tree.Size();

		assert((int)sz == lot_size*i);
		i++;
	}
 
	tree.Sync();

	struct perfmetric m;
	
	m.avg_build_ops = 100.0*((double)datapoint_t<VectorKeyObject,PL>::n_build_ops/(double)sz);
	m.avg_build_time = total.count()/(double)sz;
	
	cout << "(" << index << ") build tree: " << setw(10) << setprecision(6) << m.avg_build_ops << "% opers "
		 << setw(10) << setprecision(6) << m.avg_build_time << " microsecs ";

	double centers[n_clusters][16];
	for (int i=0;i < n_clusters;i++){
		generate_point(centers[i]);

		vector<datapoint_t<VectorKeyObject,PL>> cluster;
		generate_cluster(cluster, centers[i], clustersize, radius);
		assert((int)cluster.size() == clustersize);
		
		tree.Add(cluster);
		sz = tree.Size();

		assert((int)sz == n_points + clustersize*(i+1));
	}

	tree.Sync();

	datapoint_t<VectorKeyObject,PL>::n_query_ops = 0;
	chrono::duration<double, milli> querytime(0);
	for (int i=0;i < n_clusters;i++){
		VectorKeyObject target(centers[i]);

		auto s = chrono::steady_clock::now();
		vector<item_t<VectorKeyObject>> results = tree.Query(target, radius);
		auto e = chrono::steady_clock::now();
		querytime += (e - s);

		int nresults = (int)results.size();
		assert(nresults >= clustersize);
	}

	m.avg_query_ops = 100.0*((double)datapoint_t<VectorKeyObject,PL>::n_query_ops/(double)n_clusters)/(double)sz;
	m.avg_query_time = (double)querytime.count()/(double)n_clusters;
	m.avg_memory = tree.MemoryUsage();
	
	cout << " query ops " << dec << setprecision(6) << m.avg_query_ops << "% opers   " 
		 << "query time: " << dec << setprecision(6) << m.avg_query_time << " millisecs" << endl;

	metrics.push_back(m);

	tree.Clear();
	assert(tree.Size() == 0);
}

void do_experiment(const int n, const int n_runs, const int n_points, const int n_clusters,
				   const int clustersize, const double radius){

	cout << "-------------------Experiment " << n << " -----------------------" << endl;
	cout << "16-dim real valued vectors in L2 metric space" << endl;
	cout << "dataset size: " << n_points << " datapoints" << endl;
	cout << "no. clusters = " << n_clusters << endl;
	cout << "cluster size = " << clustersize << endl;
	cout << "query radius = " << radius << endl;
	cout << "no. test runs averaged over: " << n_runs << endl;
	
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

	cout << "no. runs: " << metrics.size() << endl;
	cout << "avg build:  " << avg_build_ops << "% opers " << avg_build_time << " microseconds" << endl;
	cout << "avg query:  " << avg_query_ops << "% opers " << avg_query_time << " milliseconds" << endl;
	cout << "avg memory: " << avg_memory/1000000 << "MB" << endl;
	cout << "-------------------------------------------------------------------" << endl << endl;
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

	const int n_runs = 5;
	const int n_experiments = 7;
	const int n_points[n_experiments] = { 100000, 200000, 400000, 800000, 1000000, 2000000, 4000000 };
	const int n_clusters = 10;
	const int clustersize = 10;
	const double radius = 0.04;

	/** 
	for (int i=0;i < n_experiments;i++){
		do_experiment(i+1, n_runs, n_points[i], n_clusters, clustersize, radius);
	}
	**/
	
	cout << "Experiments for varying radius queries" << endl << endl;
	
	const int N = 1000000;
	const int n_rad = 5;
	const double rad[n_rad] = { 0.02, 0.04, 0.06, 0.08, 0.10 };
	for (int i=0;i < n_rad;i++){
		do_experiment(i+1, n_runs, N, n_clusters, clustersize, rad[i]);
	}
	
	cout << endl << "Done." << endl;
	return 0;
}
