#include <iostream>
#include <cassert>
#include <typeinfo>
#include "mvpnode.hpp"

using namespace std;

static long long id = 1;

vector<DataPoint*> CreatePoints(const int n){
	vector<DataPoint*> results;
	for (int i=0;i<n;i++){
		DataPoint *dp = new DataPoint();
		dp->id = id++;
		dp->value = rand()%100;
		results.push_back(dp);
	}

	return results;
}

void simple_test(){
	cout << "Simple Test" << endl;

	const int n_points = 25;

	srand(1);
	vector<DataPoint*> points = CreatePoints(n_points);
	assert(points.size() == n_points);


	cout << "  Create leaf node with " << n_points << " points" << endl;
	MVPLeaf *leaf = new MVPLeaf();
	assert(leaf != NULL);

	map<int, vector<DataPoint*>*> childpoints;
	MVPNode *node = leaf->AddDataPoints(points, childpoints, 0, 0);
	assert(leaf == node);
	assert(childpoints.size() == 0);
	
	int n = node->GetCount();
	assert(n  == n_points - MVP_PATHLENGTH);
	
	DataPoint target;
	target.value = 50;
	double radius = 10;
	
	vector<DataPoint*> results = node->FilterDataPoints(&target, radius);
	assert(results.size() == 3);
	
	cout << "  Query returns " << results.size() << " results." << endl;
	for (DataPoint *dp : results){
		cout << "    id = " << dp->id << " value = " << dp->value << " "<< endl;
	}

	vector<DataPoint*> vps = node->GetVantagePoints();
	assert(vps.size() == MVP_PATHLENGTH);
	
	vector<DataPoint*> pts = node->GetDataPoints();
	assert(pts.size() == n_points - MVP_PATHLENGTH);
	

	const int n_points2 = 80;
	vector<DataPoint*> points2 = CreatePoints(n_points2);
	assert(points2.size() == n_points2);

	cout << " Add " << n_points2 << " points to leaf ==> internal" << endl;

	MVPNode *node2 = node->AddDataPoints(points2, childpoints, 0, 0);
	assert(node2 != node);

	assert(typeid(*node2).hash_code() == typeid(MVPInternal).hash_code());
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

	vector<DataPoint*> vps2 = node2->GetVantagePoints();
	assert(vps2.size() == MVP_LEVELSPERNODE);

	for (DataPoint *dp : vps2) delete dp;
	
	vector<DataPoint*> pts2 = node2->GetDataPoints();
	assert(pts2.size() == 0);

	delete node2; // delete internal node
}


int main(int argc, char **argv){

	simple_test();


	return 0;
}
