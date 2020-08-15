#include <iostream>
#include <vector>
#include <cassert>
#include "mvptree.hpp"

using namespace std;

static long long g_id = 1;

DataPoint* CreatePoint(){
	DataPoint *dp = new DataPoint();
	dp->id = g_id++;
	dp->value = rand()%100;
	return dp;
}

vector<DataPoint*> CreatePoints(const int n){
	vector<DataPoint*> results;
	for (int i=0;i<n;i++){
		DataPoint *dp = new DataPoint();
		dp->id = g_id++;
		dp->value = rand()%100;
		results.push_back(dp);
	}
	return results;
}

void simple_test(){
	const int n_points = 100;

	vector<DataPoint*> points = CreatePoints(n_points);
	assert(points.size() == n_points);

	cout << "Create MVPTree with " << n_points << " points" << endl;

	MVPTree tree;

	tree.Add(points);

	int n = tree.Size();
	assert(n == n_points);
	
	double radius = 5;
	DataPoint target;
	target.value = 50;
	vector<DataPoint*> results = tree.Query(target, radius);
	assert(results.size() == 6);

	cout << "Query for target value " << target.value << " and radius " << radius << endl;
	cout << "Found " << results.size() << " items" << endl;
	for (DataPoint *dp : results){
		cout << "id = " << dp->id << " value = " << dp->value << endl;
	}

	
	cout << "Delete item" << endl;
	tree.Delete(80);

	vector<DataPoint*> results2 = tree.Query(target, radius);
	assert(results2.size() == 5);

	cout << "Repeat Query" << endl;
	cout << "Found " << results2.size() << " items" << endl;

	size_t n_bytes = tree.MemoryUsage();
	cout << n_bytes << " bytes used" << endl;
	assert(n_bytes == 13712);
	
	cout << "Clear tree" << endl;
	tree.Clear();
	
	n = tree.Size();
	assert(n == 0);
	cout << "Number values in tree " << n << endl;
}

void test(){
	
	MVPTree tree;
	const int N = 250;
	const int n_reps = 16;

	cout << "Create MVPTree with " << n_reps << "x" << N << " points" << endl;
	size_t n_bytes;
	int n_internal, n_leaf;
	for (int n=0;n < n_reps;n++){
		for (int i=0;i<N;i++){
			DataPoint *dp = CreatePoint();
			tree.Add(dp);
		}
		tree.Sync();
		tree.CountNodes(n_internal, n_leaf);
		n_bytes = tree.MemoryUsage();
	}

	cout << "internal nodes = " << n_internal << " leaf nodes = " << n_leaf << endl;
	cout << n_bytes << " bytes used" << endl;
	assert(n_bytes == 345416);
	assert(n_internal == 233);
	assert(n_leaf == 100);
	
	int n = tree.Size();
	cout << "size = " << n << endl;
	assert(n == 4000);
	

	const DataPoint target(0, 50);
	double radius = 10;

	const vector<DataPoint*> results = tree.Query(target, radius);
	cout << "found " << results.size() << " results" << endl;
	assert(results.size() == 873);


	const DataPoint *pt = tree.Lookup(100);
	cout << "lookup point" << pt->id << " " << pt->value << endl;
	
	const DataPoint *pt2 = tree.Lookup(501);
	cout << "lookup point" << pt2->id << " " << pt2->value << endl;


	cout << "Delete " << pt->id << endl;
	tree.Delete(100);
	assert(tree.Size() == 3999);

	const DataPoint *pt3 = tree.Lookup(100);
	assert(pt3 == NULL);
	
	cout << "Delete " << pt2->id << endl;
	tree.Delete(501);
	assert(tree.Size() == 3998);

	const DataPoint *pt4 = tree.Lookup(501);
	assert(pt4 == NULL);

	cout << "Clear tree" << endl;
	tree.Clear();

	n = tree.Size();
    assert(n == 0);
	cout << "tree size = " << n << endl;

	n_bytes = tree.MemoryUsage();
	cout << "bytes = " << n_bytes << endl;
	assert(n_bytes == 2248);

	tree.CountNodes(n_internal, n_leaf);
	assert(n_internal == 0);
	assert(n_leaf == 0);
}

int main(int argc, char **argv){
	srand(1);
	
	cout << endl << "run simple test: " << endl;
	simple_test();

	g_id = 1;
	cout << endl << "run test: " << endl;
	test();

	cout << "Done." << endl;
	
	return 0;
}
