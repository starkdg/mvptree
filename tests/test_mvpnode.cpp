#include <iostream>
#include <iomanip>
#include <random>
#include <cstdint>
#include <cassert>
#include <vector>
#include <typeinfo>
#include "mvptree/mvpnode.hpp"
#include "mvptree/datapoint.hpp"

using namespace std;

static long long id = 1;
static random_device rd;
static mt19937 gen(rd());
static uniform_int_distribution<uint64_t> distrib(0);


void generate_points(vector<DataPoint*> &points, const int n){
	for (int i=0;i<n;i++){
		DataPoint *dp = new DataPoint(id++,  distrib(gen));
		points.push_back(dp);
	}
}

void test(){
	const int n_points = 25;

	vector<DataPoint*> points;
	generate_points(points, n_points);
	assert(points.size() == n_points);


	uint64_t target_value = points[0]->GetValue();

	cout << "  Create leaf node with " << n_points << " points" << endl;
	MVPLeaf<2,8,30,2,4,2> *leaf = new MVPLeaf<>();
	assert(leaf != NULL);

	map<int, vector<DataPoint*>*> childpoints;
	MVPNode<> *node = leaf->AddDataPoints(points, childpoints, 0, 0);
	assert(leaf == node);
	assert(childpoints.size() == 0);
	
	int n = node->GetCount();
	assert(n  == n_points);
	
	double radius = 5;
	
	vector<DataPoint*> results = node->FilterDataPoints(target_value, radius);
	assert(results.size() == 1);
	
	cout << "  Query returns " << results.size() << " results." << endl;
	for (DataPoint *dp : results){
		cout << "    id = " << dp->GetId() << " value = " << dp->GetValue() << endl;
	}

	vector<uint64_t> vpoints = node->GetVantagePoints();
	cout << "no. vantage points " << vpoints.size() << endl;
	assert(vpoints.size() == 8);
	
	vector<DataPoint*> pts = node->GetDataPoints();
	assert(pts.size() == n_points);
	

	const int n_points2 = 80;
	vector<DataPoint*> points2;
	generate_points(points2, n_points2);
	assert(points2.size() == n_points2);

	cout << " Add " << n_points2 << " points to leaf ==> internal" << endl;

	MVPNode<> *node2 = node->AddDataPoints(points2, childpoints, 0, 0);
	assert(node2 != node);

	assert(typeid(*node2).hash_code() == typeid(MVPInternal<>).hash_code());
	assert(childpoints.size() == 4); 
	
	delete node;

	cout << " Points sorted into lists:" << endl;
	for (auto iter=childpoints.begin();iter!=childpoints.end();iter++){
		int child_index = iter->first;
		vector<DataPoint*> *list = iter->second;
		cout << "    child[" << child_index << "] " << list->size() << " points ";
		cout << endl;
		for (DataPoint *dp : *list) delete dp;
		delete list;
	}

	vector<uint64_t> vps2 = node2->GetVantagePoints();
	assert(vps2.size() == 2);

	vector<DataPoint*> pts2 = node2->GetDataPoints();
	assert(pts2.size() == 0);

	delete node2; // delete internal node
}


int main(int argc, char **argv){

	test();


	return 0;
}
