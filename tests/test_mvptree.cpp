#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include <random>
#include "mvptree/mvptree.hpp"

using namespace std;

static long long m_id = 1;
static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);

const int BF = 2;
const int PL = 8;
const int LC = 30;
const int LPN = 2;
const int FO = 4;
const int NS = 2;


int generate_data(vector<DataPoint<PL>*> &points, int N){
	for (int i=0;i<N;i++){
		uint64_t value = m_distrib(m_gen);
		DataPoint<PL> *pnt = new H64DataPoint<PL>(m_id++, value);
		points.push_back(pnt);
	}
	return 0;
}

int generate_cluster(vector<DataPoint<PL>*> &points, uint64_t center, int N, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;

	DataPoint<PL> *mid = new H64DataPoint<PL>(m_id++, center);
	points.push_back(mid);
	for (int i=0;i < N-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		DataPoint<PL> *pnt = new H64DataPoint<PL>(m_id++, val);
		points.push_back(pnt);
	}
	return N;
}

void simple_test(){
	const int n_points = 20;

	cout << "generate " << dec << n_points << " points" << endl;
	vector<DataPoint<PL>*> points;
	generate_data(points, n_points);

	assert(points.size() == n_points);

	cout << "build tree" << endl;
	MVPTree<BF,PL,LC,LPN,FO,NS> tree;
	tree.Add(points);
	tree.Sync();
	
	int n = tree.Size();
	cout << "tree size is " << dec << n << endl;
	assert(n == n_points);

	uint64_t center = m_distrib(m_gen);
	double radius = 5;
	int cluster_size = 5;
	vector<DataPoint<PL>*> cluster;

	cout << "generate cluster of " << dec << cluster_size << " points" << endl;
	generate_cluster(cluster, center, cluster_size, radius);
	assert((int)cluster.size() == cluster_size);
	tree.Add(cluster);
	tree.Sync();
	
	n = tree.Size();
	cout << "tree size is " << dec << n << endl;
	assert(n == n_points + cluster_size);
	
	
	H64DataPoint<PL> target(0, center);
	vector<DataPoint<PL>*> results = tree.Query(target, radius);
	assert((int)results.size() == cluster_size);

	cout << "Query for target value " << target.GetValue() << " and radius " << radius << endl;
	cout << "Found " << results.size() << " items" << endl;
	assert((int)results.size() == cluster_size);
	for (DataPoint<PL> *dp : results){
		H64DataPoint<PL> *pnt = (H64DataPoint<PL>*)dp;
		cout << "==> id = " << pnt->GetId() << " value = " << hex << pnt->GetValue() << endl;
	}

	
	cout << "Delete item" << endl;
	tree.Delete(5);

	n = tree.Size();
	cout << dec << "tree size: " << n << endl;
	assert(n == n_points + cluster_size - 1);
	
	size_t n_bytes = tree.MemoryUsage();
	cout << dec << n_bytes << " bytes used" << endl;
	assert(n_bytes == 2872);


	const int n_points2 = 30;

	cout << "generate " << dec << n_points2 << " points" << endl;
	vector<DataPoint<PL>*> points2;
	generate_data(points2, n_points2);

	cout << "Add " << dec << n_points2 << " to tree" << endl;
	tree.Add(points2);
	n = tree.Size();
	cout << "tree size " << dec << n << endl;
	assert(n == n_points + cluster_size + n_points2 - 1);


	int internal_nodes = 0, leaf_nodes = 0;
	tree.CountNodes(internal_nodes, leaf_nodes);
	cout << "internal: " << internal_nodes << " leaf: " << leaf_nodes << endl;
	assert(internal_nodes == 1 && leaf_nodes == 4);
	
	const H64DataPoint<PL> *retnode = (H64DataPoint<PL>*)tree.Lookup(10);
	assert(retnode != NULL);
	cout << "Look up id=10 => " << retnode->GetId() << " " << hex << retnode->GetValue() << endl;
	
	cout << "Clear tree" << endl;
	tree.Clear();
	
	n = tree.Size();
	cout << "tree size " << n << endl;
	assert(n == 0);
}

void test(){
	const int N = 100;
	const int n_iters = 5;

	MVPTree<BF,PL,LC,LPN,FO,NS> mvptree;
	for (int i=0;i < n_iters;i++){
		vector<DataPoint<PL>*> points;
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
	const int radius = 10;
	uint64_t centers[n_clusters];
	for (int i=0;i < n_clusters;i++){
		cout << "Add cluster of " << cluster_size << " points to tree" << endl;

		centers[i] = m_distrib(m_gen);

		vector<DataPoint<PL>*> cluster;
		generate_cluster(cluster, centers[i], cluster_size, radius);
		mvptree.Add(cluster);
	}
		
	mvptree.Sync();
	sz = mvptree.Size();
	cout << "tree size: " << sz << endl;
	assert(sz == N*n_iters + cluster_size*n_clusters);

	for (int i=0;i < n_clusters;i++){
		H64DataPoint<PL> target(0, centers[i]);
		cout << "Query " << target.GetId() << " " << hex << target.GetValue() << endl;
		vector<DataPoint<PL>*> results = mvptree.Query(target, radius);
		cout << "==> Found " << dec << results.size() << " points" << endl;
		assert(results.size() >= n_clusters);
	}

	int n_internal, n_leaf;
	mvptree.CountNodes(n_internal, n_leaf);
	cout << "internal nodes: " << n_internal << endl;
	cout << "leaf nodes: " << n_leaf << endl;

	
	const DataPoint<PL> *dp = mvptree.Lookup(100);
	cout << "Lookup 100 => " << dp->GetId() << endl;
	assert(dp != NULL);
	
	const DataPoint<PL> *dp2 = mvptree.Lookup(155);
	cout << "Lookup 155 => " << dp2->GetId() << endl;
	assert(dp2 != NULL);

	cout << "Delete 100" << endl;
	mvptree.Delete(100);

	cout << "Delete 155" << endl;
	mvptree.Delete(155);

	sz = mvptree.Size();
	cout << "tree size: " << sz << endl;
	assert(sz == N*n_iters + cluster_size*n_clusters - 2);

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
