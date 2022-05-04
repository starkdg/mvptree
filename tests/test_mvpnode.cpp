#include <iostream>
#include <cassert>
#include <cstdint>
#include <typeinfo>
#include <random>
#include "mvptree/mvpnode.hpp"
#include "mvptree/datapoint.hpp"
#include "mvptree/key.hpp"

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


int generate_data(vector<datapoint_t<H64KeyObject,PL>> &points, int N){
	for (int i=0;i<N;i++){
		uint64_t value = m_distrib(m_gen);
		datapoint_t<H64KeyObject,PL> pnt(m_id++, H64KeyObject(value));
		points.push_back(pnt);
	}
	return 0;
}

int generate_cluster(vector<datapoint_t<H64KeyObject,PL>> &points, uint64_t center, int N, int max_radius){
	static uniform_int_distribution<int> radius_distr(1, max_radius);
	static uniform_int_distribution<int> bitindex_distr(0, 63);
		
	uint64_t mask = 0x01;

	datapoint_t<H64KeyObject,PL> mid(m_id++, H64KeyObject(center));
	points.push_back(mid);
	for (int i=0;i < N-1;i++){
		uint64_t val = center;
		int dist = radius_distr(m_gen);
		for (int j=0;j < dist;j++){
			val ^= (mask << bitindex_distr(m_gen));
		}
		datapoint_t<H64KeyObject,PL> pnt(m_id++, H64KeyObject(val));
		points.push_back(pnt);
	}
	return N;
}

void test_mvpnode(){
	const int n_points = 20;

	vector<datapoint_t<H64KeyObject,PL>> points;
	generate_data(points, n_points);
	assert(points.size() == n_points);

	cout << "Create Leaf node for points" << endl;
	MVPLeaf<H64KeyObject,BF,PL,LC,LPN,FO,NS> *leaf = new MVPLeaf<H64KeyObject,BF,PL,LC,LPN,FO,NS>();
	assert(leaf != NULL);

	map<int, vector<datapoint_t<H64KeyObject,PL>>*>  childpoints;
	MVPNode<H64KeyObject,BF,PL,LC,LPN,FO,NS> *node = leaf->AddDataPoints(points, childpoints, 0, 0);
	assert(leaf == node);
	assert(childpoints.size() == 0);
	
	int n = node->GetCount();
	cout << "leaf contains " << n << " points" << endl;
	assert(n  == n_points - 2);


	const int cluster_size = 5;
	const int radius = 10;
	vector<datapoint_t<H64KeyObject,PL>> cluster;

	uint64_t center = m_distrib(m_gen);
	generate_cluster(cluster, center, cluster_size, radius);

	
	MVPNode<H64KeyObject,BF,PL,LC,LPN,FO,NS> *node2 = node->AddDataPoints(cluster, childpoints, 0, 0);
	assert(node2 == node);


	H64KeyObject target(center);
	
	vector<double> path;
	vector<item_t<H64KeyObject>> results;
	node->FilterDataPoints(target, path, radius, 0, false, results);
	cout << "Query target: " << hex << center << endl;
	cout << "Found: " << dec << results.size() << endl;
	assert(results.size() == cluster_size);
	
	vector<vp_t<H64KeyObject>> vps = node->GetVantagePoints();
	cout << "vantage points: " << vps.size() << endl;
	assert(vps.size() == 2);
	
	vector<datapoint_t<H64KeyObject,PL>> pts = node->GetDataPoints();
	cout << "no. points: " << pts.size() << endl;
	assert(pts.size() == n_points + cluster_size - 2);

	
	const int n_points2 = 80;
	vector<datapoint_t<H64KeyObject,PL>> points2;
	generate_data(points2, n_points2);
	assert(points2.size() == n_points2);

	cout << " Add " << n_points2 << " points to leaf ==> internal" << endl;

	MVPNode<H64KeyObject,BF,PL,LC,LPN,FO,NS> *internal = node->AddDataPoints(points2, childpoints, 0, 0);
	assert(internal != node);

	assert(typeid(*internal).hash_code() == typeid(MVPInternal<H64KeyObject,BF,PL,LC,LPN,FO,NS>).hash_code());
	cout << "child nodes: " << childpoints.size() << endl;
	assert(childpoints.size() == 4);

	
	delete node;

	for (auto iter=childpoints.begin();iter!=childpoints.end();iter++){
		int child_index = iter->first;
		vector<datapoint_t<H64KeyObject,PL>> *list = iter->second;
		cout << "    child[" << child_index << "] " << list->size() << " points ";
		cout << endl;
		delete list;
	}

	const vector<vp_t<H64KeyObject>> vps2 = internal->GetVantagePoints();
	assert(vps2.size() == 2);

	const vector<datapoint_t<H64KeyObject,PL>> pts2 = internal->GetDataPoints();
	assert(pts2.size() == 0);

	delete internal;
}


int main(int argc, char **argv){

	cout << "Run test" << endl;
	test_mvpnode();


	return 0;
}
