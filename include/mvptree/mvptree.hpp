#ifndef _MVPTREE_H
#define _MVPTREE_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <queue>
#include "mvptree/mvpnode.hpp"
#include "mvptree/datapoint.hpp"

using namespace std;

template<int BF, int PL, int LC, int LPN, int FO, int NS>
class MVPNode;

template<int BF, int PL, int LC, int LPN, int FO, int NS>
class MVPInternal;

template<int BF, int PL, int LC, int LPN, int FO, int NS>
class MVPLeaf;

/**
 * template arguments
 *  BF - branchfactor, no. of units to split space at each level of tree
 *  PL - pathlength, no. vantage point distances to maintain for each datapoint
 *  LC - leafcapacity, no. points in each leaf node
 *  LPN - levelspernode, no. levels per internal node in which the space is successively split into BF
 *  FO - fanout, no childnodes per internal node BF^LPN
 *  NS - numsplits, max no. split points for last level of internal node BF^(LPN-1)
 **/

template<int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
class MVPTree {
private:
	vector<DataPoint<PL>*> m_arrivals;
	
	map<long long, DataPoint<PL>*> m_ids;
	
	MVPNode<BF,PL,LC,LPN,FO,NS> *m_top;

	int n_internal, n_leaf, n_sync;

	void LinkNodes(map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &nodes,
				   map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes)const;
	void ExpandNode(MVPNode<BF,PL,LC,LPN,FO,NS> *node,
					map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes, const int index)const;
	MVPNode<BF,PL,LC,LPN,FO,NS>* ProcessNode(const int level,
											 const int index,
											 MVPNode<BF,PL,LC,LPN,FO,NS> *node, vector<DataPoint<PL>*> &points,
											 map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
											 map<int, vector<DataPoint<PL>*>*> &childpoints);
public:
	MVPTree():m_top(NULL),n_internal(0),n_leaf(0),n_sync(100){};

	const DataPoint<PL>* Lookup(const long long id);
	
	void Add(DataPoint<PL> *dp);
	
	void Add(vector<DataPoint<PL>*> &points);

	void Sync();
	
	void Delete(const long long id);

	const int Size()const;

	void CountNodes(int &n_internal, int &n_leaf)const;
	
	void Clear();

	const vector<DataPoint<PL>*> Query(const DataPoint<PL> &target, const double radius) const;

	size_t MemoryUsage()const;

	void SetSync(int n);
};


/** -----------------------------------
 *
 *     MVPTree implementation
 *
 **/

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::LinkNodes(map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &nodes,
											map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes)const{
	for (auto iter=nodes.begin();iter!=nodes.end();iter++){
		int i = iter->first;
		MVPNode<BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;
		if (mvpnode != NULL){
			for (int j=0;j < FO;j++){
				MVPNode<BF,PL,LC,LPN,FO,NS> *child = childnodes[i*FO+j];
				if (child != NULL) mvpnode->SetChildNode(j, child);
			}
		}
	}
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::ExpandNode(MVPNode<BF,PL,LC,LPN,FO,NS> *node,
											 map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
											 const int index)const{
	if (node != NULL){
		for (int i=0;i < FO;i++){
			MVPNode<BF,PL,LC,LPN,FO,NS> *child = node->GetChildNode(i);
			if (child != NULL) childnodes[index*FO+i] = child;
		}
	}
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPTree<BF,PL,LC,LPN,FO,NS>::ProcessNode(const int level, const int index,
							  MVPNode<BF,PL,LC,LPN,FO,NS> *node,
							  vector<DataPoint<PL>*> &points,
							  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
							  map<int, vector<DataPoint<PL>*>*> &childpoints){
	MVPNode<BF,PL,LC,LPN,FO,NS> *retnode = node;
	if (node == NULL){ // create new node
		retnode = MVPNode<BF,PL,LC,LPN,FO,NS>::CreateNode(points, childpoints, level, index);
	} else {           // node exists
		retnode = node->AddDataPoints(points, childpoints, level, index);
	}

	if (retnode == NULL)
		throw runtime_error("failure to assign node");
	
	return retnode;
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
const DataPoint<PL>* MVPTree<BF,PL,LC,LPN,FO,NS>::Lookup(const long long id){
	auto iter = m_ids.find(id);
	if (iter != m_ids.end()){
		return iter->second;
	}
	return NULL;
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::Add(DataPoint<PL> *dp){
	if (dp != NULL){
		m_arrivals.push_back(dp);
		if ((int)m_arrivals.size() >= n_sync) Add(m_arrivals);
	}
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::Add(vector<DataPoint<PL>*> &points){
	if (points.empty()) return;

	for (DataPoint<PL>* dp : points) m_ids[dp->GetId()] = dp;

	map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> prevnodes, currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	map<int, vector<DataPoint<PL>*>*> pnts, pnts2;
	pnts[0] = &points;

	int n = 0;
	do {
		for (auto iter=pnts.begin();iter!=pnts.end();iter++){
			int index = iter->first;
			vector<DataPoint<PL>*> *list = iter->second;
			MVPNode<BF,PL,LC,LPN,FO,NS> *mvpnode = currnodes[index];
			MVPNode<BF,PL,LC,LPN,FO,NS> *newnode = ProcessNode(n, index, mvpnode, *list, childnodes, pnts2);
			if (newnode != mvpnode){
				if (newnode != NULL && typeid(*newnode).hash_code() == typeid(MVPInternal<BF,PL,LC,LPN,FO,NS>).hash_code()){
					n_internal++;
				} else if (typeid(*newnode).hash_code() == typeid(MVPLeaf<BF,PL,LC,LPN,FO,NS>).hash_code()){
					n_leaf++;
				}
				if (mvpnode != NULL){
					if (typeid(*mvpnode).hash_code() == typeid(MVPLeaf<BF,PL,LC,LPN,FO,NS>).hash_code()){
						n_leaf--;
					}
					delete mvpnode;
				}
				if (n == 0) m_top = newnode;
			}
			currnodes[index] = newnode;
			ExpandNode(newnode, childnodes, index);

			if (n > 0) delete list;
		}

		if (!prevnodes.empty()) {
			LinkNodes(prevnodes, currnodes);
		}
		prevnodes = move(currnodes);
		currnodes = move(childnodes);
		pnts = move(pnts2);
		childnodes.clear();
		n += LPN;
	} while (!pnts.empty());
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::Sync(){
	if (m_arrivals.size() > 0) {
		Add(m_arrivals);
	}
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::Delete(const long long id){
	auto iter = m_ids.find(id);
	if (iter != m_ids.end()) iter->second->Deactivate();
	m_ids.erase(id);
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
const int MVPTree<BF,PL,LC,LPN,FO,NS>::Size()const{
	return m_ids.size();
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::CountNodes(int &n_internal, int &n_leaf)const{
	n_internal = 0;
	n_leaf = 0;

	queue<MVPNode<BF,PL,LC,LPN,FO,NS>*> nodes;
	if (m_top != NULL) nodes.push(m_top);

	while (!nodes.empty()){
		MVPNode<BF,PL,LC,LPN,FO,NS> *curr_node = nodes.front();

		if (typeid(*curr_node).hash_code() == typeid(MVPInternal<BF,PL,LC,LPN,FO,NS>).hash_code()){
			n_internal++;
		} else {
			n_leaf++;
		}
	
		for (int i=0;i < FO;i++){
			MVPNode<BF,PL,LC,LPN,FO,NS> *child = curr_node->GetChildNode(i);
			if (child != NULL) nodes.push(child);
		}
		
		nodes.pop();
	}
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::Clear(){
	map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n = 0, n_internal = 0, n_leaf = 0;
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
			int index = iter->first;
			MVPNode<BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;

			if (mvpnode != NULL){
				vector<DataPoint<PL>*> pts = mvpnode->PurgeDataPoints();
				for (DataPoint<PL> *dp : pts){
					delete dp;
				}

				ExpandNode(mvpnode, childnodes, index);
				
				if (typeid(*mvpnode).hash_code() == typeid(MVPInternal<BF,PL,LC,LPN,FO,NS>).hash_code()){
					n_internal++;
				} else {
					n_leaf++;
				}

				delete mvpnode;
			}
		}
		currnodes = move(childnodes);
		n += LPN;

	} while (!currnodes.empty());
	m_top = NULL;
	m_ids.clear();
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
const vector<DataPoint<PL>*> MVPTree<BF,PL,LC,LPN,FO,NS>::Query(const DataPoint<PL> &target, const double radius) const{
	vector<DataPoint<PL>*> results;
	map<int, vector<double>*> tdistances, tdistances2;
	map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;
	
	int n = 0;
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		
		for (auto const &iter : currnodes){
			int node_index = iter.first;
			MVPNode<BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
			vector<double> *tpath = tdistances[node_index];
			if (mvpnode)
				mvpnode->TraverseNode(target, radius, childnodes,  tpath, tdistances2, node_index, n, results);
		}
		currnodes = move(childnodes);
		tdistances = move(tdistances2);
		n += LPN;
	} while (!currnodes.empty());

	return results;
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
size_t MVPTree<BF,PL,LC,LPN,FO,NS>::MemoryUsage()const{
	map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n_points = m_ids.size();
	int n_internal=0, n_leaf=0;
	do {
		for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
			int index = iter->first;
			MVPNode<BF,PL,LC,LPN,FO,NS> *node = iter->second;

			if (typeid(*node).hash_code() == typeid(MVPInternal<BF,PL,LC,LPN,FO,NS>).hash_code())
				n_internal++;
			else
				n_leaf++;
			ExpandNode(node, childnodes, index);
		}
		currnodes = move(childnodes);
	} while (!currnodes.empty());
	
	return  n_points*sizeof(DataPoint<PL>) + n_internal*sizeof(MVPInternal<BF,PL,LC,LPN,FO,NS>)
		+ n_leaf*sizeof(MVPLeaf<BF,PL,LC,LPN,FO,NS>) + sizeof(MVPLeaf<BF,PL,LC,LPN,FO,NS>);
}

template<int BF,int PL, int LC, int LPN, int FO, int NS>
void MVPTree<BF,PL,LC,LPN,FO,NS>::SetSync(int n){
	n_sync = n;
}


#endif /* _MVPTREE_H */
