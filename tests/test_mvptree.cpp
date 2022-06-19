#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include <random>
#include "mvptree/mvptree.hpp"
#include "mvptree/key.hpp"

using namespace std;
using namespace mvp;

static long long m_id = 1;
static long long g_id = 1000;

static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);

const int BF = 2;
const int PL = 8;
const int LC = 30;
const int LPN = 2;
const int FO = 4;
const int NS = 2;



/**
 * generate N datapoints
 **/
int generate_data(vector<datapoint_t<H64KeyObject,PL>> &points, int N){
	for (int i=0;i<N;i++){
		uint64_t value = m_distrib(m_gen);
		points.push_back({ m_id++, H64KeyObject(value) }); 
	}
	return 0;
}

/**
 * generate a cluster of N datapoints with center point and cluster size of max_radius
 **/
int generate_cluster(vector<datapoint_t<H64KeyObject,PL>> &points, uint64_t center, int N, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;
	points.push_back({g_id++, H64KeyObject(center) });
	for (int i=0;i < N-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		points.push_back({ g_id++, H64KeyObject(val)} );
	}
	return N;
}

void simple_test(){
	const int n_points = 20;

	cout << "generate " << dec << n_points << " points" << endl;
	vector<datapoint_t<H64KeyObject,PL>> points;
	generate_data(points, n_points);

	assert(points.size() == n_points);

	cout << "build tree" << endl;
	MVPTree<H64KeyObject, BF,PL,LC,LPN,FO,NS> tree;
	tree.Add(points);
	tree.Sync();
	
	int n = tree.Size();
	cout << "tree size is " << dec << n << endl;
	assert(n == n_points);

	uint64_t center = m_distrib(m_gen);
	double radius = 5;
	int cluster_size = 5;
	vector<datapoint_t<H64KeyObject,PL>> cluster;

	cout << "generate cluster of " << dec << cluster_size << " points" << endl;
	generate_cluster(cluster, center, cluster_size, radius);
	assert((int)cluster.size() == cluster_size);
	tree.Add(cluster);
	tree.Sync();
	
	n = tree.Size();
	cout << "tree size is " << dec << n << endl;
	assert(n == n_points + cluster_size);
	
	H64KeyObject target(center);
	vector<item_t<H64KeyObject>> results = tree.Query(target, radius);
	assert((int)results.size() == cluster_size);

	cout << "Query for target " << center << " and radius " << radius << endl;
	cout << "Found " << results.size() << " items" << endl;
	assert((int)results.size() == cluster_size);
	for (item_t<H64KeyObject> &item : results){
		cout << "==> id = " << item.id << " value = " << hex << item.key.key << endl;
	}

	
	n = tree.Size();
	cout << dec << "tree size: " << n << endl;
	assert(n == n_points + cluster_size);
	
	size_t n_bytes = tree.MemoryUsage();
	cout << dec << n_bytes << " bytes used" << endl;
	//assert(n_bytes == 4208);
	
	const int n_points2 = 30;

	cout << "generate " << dec << n_points2 << " points" << endl;
	vector<datapoint_t<H64KeyObject,PL>> points2;
	generate_data(points2, n_points2);

	cout << "Add " << dec << n_points2 << " to tree" << endl;
	tree.Add(points2);
	n = tree.Size();
	cout << "tree size " << dec << n << endl;
	assert(n == n_points + cluster_size + n_points2);


	cout << "Clear tree" << endl;
	tree.Clear();
	
	n = tree.Size();
	cout << "tree size " << n << endl;
	assert(n == 0);
}

void test(){
	const int N = 100;
	const int n_iters = 5;

	MVPTree<H64KeyObject,BF,PL,LC,LPN,FO,NS> mvptree;
	for (int i=0;i < n_iters;i++){
		vector<datapoint_t<H64KeyObject,PL>> points;
		cout << "Add " << dec << N << " points" << endl;
		generate_data(points, N);
		mvptree.Add(points);
		int count = mvptree.Size();
		assert(count == N*(i+1));
		
	}
	mvptree.Sync();
	
	int sz = mvptree.Size();
	cout << "tree size: " << sz << endl;
	assert(sz == N*n_iters);


	const int cluster_size = 10;
	const int n_clusters = 5;
	const double radius = 10;
	uint64_t centers[n_clusters];
	for (int i=0;i < n_clusters;i++){
		cout << "Add cluster of " << cluster_size << " points to tree" << endl;

		centers[i] = m_distrib(m_gen);

		vector<datapoint_t<H64KeyObject,PL>> cluster;
		generate_cluster(cluster, centers[i], cluster_size, radius);
		mvptree.Add(cluster);
	}
		
	mvptree.Sync();
	sz = mvptree.Size();
	cout << "tree size: " << sz << endl;
	assert(sz == N*n_iters + cluster_size*n_clusters);

	for (int i=0;i < n_clusters;i++){
		H64KeyObject  target(centers[i]);
		cout << "Query " << hex << centers[i] << endl;
		vector<item_t<H64KeyObject>> results = mvptree.Query(target, radius);
		cout << "==> Found " << dec << results.size() << " points" << endl;
		assert(results.size() >= n_clusters);
	}

	cout << "Delete id == 1000" << endl;
	int ndels = mvptree.DeletePoint(H64KeyObject(centers[0]));
	cout << "Deleted " << ndels << " points" << endl;
	assert(ndels == 1);

	sz = mvptree.Size();
	cout << "tree size: " << sz << endl;
	assert(sz == N*n_iters + cluster_size*n_clusters - 1);

	cout << "Clear tree" << endl;
	mvptree.Clear();
	sz = mvptree.Size();
	assert(sz == 0);
	
}

int main(int argc, char **argv){
	
	cout << endl << "run simple test: " << endl;
	simple_test();


	cout << endl << "run test:" << endl;
	test();
	
	cout << "Done." << endl;
	
	return 0;
}
