#include <iostream>
#include <cassert>
#include <cstdint>
#include <typeinfo>
#include <random>
#include "mvptree/mvpnode.hpp"
#include "mvptree/datapoint.hpp"

using namespace std;

static long long m_id = 1;
static random_device m_rd;
static mt19937_64 m_gen(m_rd());
static uniform_int_distribution<uint64_t> m_distrib(0);


int generate_data(vector<DataPoint<>*> &points, int N){
	for (int i=0;i<N;i++){
		uint64_t value = m_distrib(m_gen);
		DataPoint<> *pnt = new H64DataPoint<>(m_id++, value);
		points.push_back(pnt);
	}
	return 0;
}

int generate_cluster(vector<DataPoint<>*> &points, uint64_t center, int N, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;

	DataPoint<> *mid = new H64DataPoint<>(m_id++, center);
	points.push_back(mid);
	for (int i=0;i < N-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		DataPoint<> *pnt = new H64DataPoint<>(m_id++, val);
		points.push_back(pnt);
	}
	return N;
}

void test_mvpnode(){
	const int n_points = 20;

	vector<DataPoint<>*> points;
	generate_data(points, n_points);
	assert(points.size() == n_points);

	cout << "Create Leaf node for points" << endl;
	MVPLeaf<> *leaf = new MVPLeaf<>();
	assert(leaf != NULL);

	map<int, vector<DataPoint<>*>*> childpoints;
	MVPNode<> *node = leaf->AddDataPoints(points, childpoints, 0, 0);
	assert(leaf == node);
	assert(childpoints.size() == 0);
	
	int n = node->GetCount();
	cout << "leaf contains " << n << " points" << endl;
	assert(n  == n_points - 2);


	const int cluster_size = 5;
	const int radius = 10;
	vector<DataPoint<>*> cluster;

	uint64_t center = m_distrib(m_gen);
	generate_cluster(cluster, center, cluster_size, radius);

	
	MVPNode<> *node2 = node->AddDataPoints(cluster, childpoints, 0, 0);
	assert(node2 == node);
	
	
	
	const H64DataPoint<> target(0, center);
	
	vector<double> path;
	vector<DataPoint<>*> results;
	node->FilterDataPoints(target, path, radius, 0, results);
	cout << "Query target: " << hex << center << endl;
	cout << "Found: " << dec << results.size() << endl;
	assert(results.size() == cluster_size);
	
	vector<DataPoint<>*> vps = node->GetVantagePoints();
	cout << "vantage points: " << vps.size() << endl;
	assert(vps.size() == 2);
	
	vector<DataPoint<>*> pts = node->GetDataPoints();
	cout << "no. points: " << pts.size() << endl;
	assert(pts.size() == n_points + cluster_size - 2);

	
	const int n_points2 = 80;
	vector<DataPoint<>*> points2;
	generate_data(points2, n_points2);
	assert(points2.size() == n_points2);

	cout << " Add " << n_points2 << " points to leaf ==> internal" << endl;

	MVPNode<> *internal = node->AddDataPoints(points2, childpoints, 0, 0);
	assert(internal != node);

	assert(typeid(*internal).hash_code() == typeid(MVPInternal<>).hash_code());
	cout << "child nodes: " << childpoints.size() << endl;
	assert(childpoints.size() == 4);

	
	delete node;

	for (auto iter=childpoints.begin();iter!=childpoints.end();iter++){
		int child_index = iter->first;
		vector<DataPoint<>*> *list = iter->second;
		cout << "    child[" << child_index << "] " << list->size() << " points ";
		cout << endl;
		for (int i = 0;i < (int)list->size();i++)
			delete list->at(i);
		delete list;
	}

	const vector<DataPoint<>*> vps2 = internal->GetVantagePoints();
	assert(vps2.size() == 2);

	for (DataPoint<> *dp: vps2)
		delete dp;
	
	const vector<DataPoint<>*> pts2 = internal->GetDataPoints();
	assert(pts2.size() == 0);

	delete internal;
}


int main(int argc, char **argv){

	cout << "Run test" << endl;
	test_mvpnode();


	return 0;
}
