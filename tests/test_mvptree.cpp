#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include <random>
#include "mvptree/mvptree.hpp"

using namespace std;

static long long g_id = 1;
static random_device rd;
static mt19937 gen(rd());
static uniform_int_distribution<uint64_t> distrib(0);


void generate_points(vector<DataPoint*> &points, const int n){
	for (int i=0;i<n;i++){
		DataPoint *dp = new DataPoint(g_id++, distrib(gen));
		points.push_back(dp);
	}
}

void test(){
	const int n_points = 100;

	vector<DataPoint*> points;
	generate_points(points, n_points);
	assert(points.size() == n_points);

	uint64_t target_value = points[0]->GetValue();
	
	MVPTree<2,8,30,2,4,2> tree;
	tree.Add(points);
	tree.Sync();

	
	int n = tree.Size();
	assert(n == n_points);
	
	double radius = 5;
	vector<DataPoint*> results = tree.Query(target_value, radius);
	assert(results.size() == 1);
	
	cout << "Query for target value " << target_value << " and radius " << radius << endl;
	cout << "Found " << results.size() << " items" << endl;
	cout << "    id = " << results[0]->GetId() << " value = " << results[0]->GetId() << endl;
	assert(results[0]->GetId() == 1);
	
	
	cout << "Delete item" << endl;
	tree.Delete(80);


	n = tree.Size();
	assert(n == n_points - 1);
	
	size_t n_bytes = tree.MemoryUsage();
	assert(n_bytes == 13712);
	
	cout << n_bytes << " bytes used" << endl;
	
	cout << "Clear tree" << endl;
	tree.Clear();
	n = tree.Size();
	cout << "no. points " << n << endl;
	assert(n == 0);
}


int main(int argc, char **argv){

	test();
	
	return 0;
}
