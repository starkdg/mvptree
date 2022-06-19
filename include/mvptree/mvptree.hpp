#ifndef _MVPTREE_H
#define _MVPTREE_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <queue>
#include <ostream>
#include "mvptree/mvpnode.hpp"
#include "mvptree/datapoint.hpp"


/**
 * template arguments
 *  T  - key type defined in  key.hpp
 *  BF - branchfactor, no. of units to split space at each level of tree
 *  PL - pathlength, no. vantage point distances to maintain for each datapoint
 *  LC - leafcapacity, no. points in each leaf node
 *  LPN - levelspernode, no. levels per internal node in which the space is successively split into BF
 *  FO - fanout, no childnodes per internal node BF^LPN
 *  NS - numsplits, max no. split points for last level of internal node BF^(LPN-1)
 **/

namespace mvp {

	template<typename T,int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
	class MVPTree {
	private:
		std::vector<datapoint_t<T,PL>> m_arrivals;
	
		MVPNode<T,BF,PL,LC,LPN,FO,NS> *m_top;

		int n_sync;

		void LinkNodes(std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &nodes,
					   std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes)const;
		void ExpandNode(MVPNode<T,BF,PL,LC,LPN,FO,NS> *node,
						std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes, const int index)const;
		MVPNode<T,BF,PL,LC,LPN,FO,NS>* ProcessNode(const int level,
												   const int index,
												   MVPNode<T,BF,PL,LC,LPN,FO,NS> *node, std::vector<datapoint_t<T,PL>> &points,
												   std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
												   std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints);
	public:
		MVPTree():m_top(NULL),n_sync(100){};

		void Add(datapoint_t<T,PL> &item);
	
		void Add(std::vector<datapoint_t<T,PL>> &items);
		
		void Sync();
	
		const size_t Size()const;

		void Clear();

		const std::vector<item_t<T>> Query(const T &target, const double radius) const;

		const int DeletePoint(const T &target);
	
		void Print(std::ostream &ostrm)const;

		size_t MemoryUsage()const;

		void SetSync(int n);
	};
}

/** -----------------------------------
 *
 *     MVPTree implementation
 *
 **/

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::LinkNodes(std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &nodes,
											std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes)const{
	for (auto iter=nodes.begin();iter!=nodes.end();iter++){
		int i = iter->first;
		MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;
		if (mvpnode != NULL){
			for (int j=0;j < FO;j++){
				MVPNode<T,BF,PL,LC,LPN,FO,NS> *child = childnodes[i*FO+j];
				if (child != NULL) mvpnode->SetChildNode(j, child);
			}
		}
	}
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::ExpandNode(MVPNode<T,BF,PL,LC,LPN,FO,NS> *node,
											 std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
											 const int index)const{
	if (node != NULL){
		for (int i=0;i < FO;i++){
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *child = node->GetChildNode(i);
			if (child != NULL) childnodes[index*FO+i] = child;
		}
	}
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::ProcessNode(const int level,
																		  const int index,
																		  MVPNode<T,BF,PL,LC,LPN,FO,NS> *node,
																		  std::vector<datapoint_t<T,PL>> &points,
																		  std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
																		  std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints){
	MVPNode<T,BF,PL,LC,LPN,FO,NS> *retnode = node;
	if (node == NULL){ // create new node
		retnode = MVPNode<T,BF,PL,LC,LPN,FO,NS>::CreateNode(points, childpoints, level, index);
	} else {           // node exists
		retnode = node->AddDataPoints(points, childpoints, level, index);
	}

	if (retnode == NULL)
		throw std::runtime_error("failure to assign node");
	
	return retnode;
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Add(datapoint_t<T,PL> &dp){
	m_arrivals.push_back({ dp.id, dp.key });
	if ((int)m_arrivals.size() >= n_sync){
		Add(m_arrivals);
	}
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Add(std::vector<datapoint_t<T,PL>> &points){
	if (points.empty()) return;

	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> prevnodes, currnodes, childnodes;
	if (m_top != NULL)
		currnodes[0] = m_top;

	std::map<int, std::vector<datapoint_t<T,PL>>*> pnts, pnts2;
	pnts[0] = &points;

	int n = 0;
	do {
		for (auto iter=pnts.begin();iter!=pnts.end();iter++){
			int index = iter->first;
			std::vector<datapoint_t<T,PL>> *list = iter->second;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = currnodes[index];
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *newnode = ProcessNode(n, index, mvpnode, *list, childnodes, pnts2);
			if (newnode != mvpnode){
				if (mvpnode != NULL){
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

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Sync(){
	if (m_arrivals.size() > 0) {
		Add(m_arrivals);
	}
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
const size_t mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Size()const{
	size_t count = 0;
	
   	std::queue<MVPNode<T,BF,PL,LC,LPN,FO,NS>*> nodes;
	if (m_top != NULL) nodes.push(m_top);

	while (!nodes.empty()){
		MVPNode<T,BF,PL,LC,LPN,FO,NS> *current = nodes.front();

		count += current->Size();
	
		for (int i=0;i < FO;i++){
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *child = current->GetChildNode(i);
			if (child != NULL) nodes.push(child);
		}
		
		nodes.pop();
	}
	return count;
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Clear(){
	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n = 0;
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
			int index = iter->first;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;

			if (mvpnode != NULL){
				ExpandNode(mvpnode, childnodes, index);
				delete mvpnode;
			}
		}
		currnodes = move(childnodes);
		n += LPN;

	} while (!currnodes.empty());
	m_top = NULL;
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
const std::vector<mvp::item_t<T>> mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Query(const T &target, const double radius) const{
	std::vector<item_t<T>> results;
	std::map<int, std::vector<double>*> tdistances, tdistances2;
	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;
	
	int n = 0;
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		
		for (auto const &iter : currnodes){
			int node_index = iter.first;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
			std::vector<double> *tpath = tdistances[node_index];
			if (mvpnode)
				mvpnode->TraverseNode(target, radius, childnodes,  tpath, tdistances2, node_index, n, false, results);
		}
		currnodes = move(childnodes);
		tdistances = move(tdistances2);
		n += LPN;
	} while (!currnodes.empty());

	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::DeletePoint(const T &target){
	std::vector<item_t<T>> results;

	std::map<int, std::vector<double>*> tdistances, tdistances2;
	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n = 0;
	while (!currnodes.empty()) {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		for (auto const &iter : currnodes){
			int node_index = iter.first;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
			std::vector<double> *tpath = tdistances[node_index];
			if (mvpnode){
				mvpnode->TraverseNode(target, 0, childnodes, tpath, tdistances2, node_index, n, true, results);
			}
		}
		currnodes = move(childnodes);
		tdistances = move(tdistances2);
		n += LPN;
	};

	return (int)results.size();
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::Print(std::ostream &ostrm)const {
	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n=0;
	do {

		ostrm << std::dec << "level = " << n << " (" << currnodes.size() << " nodes" << std::endl;;
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+LPN);
		for (auto const &iter : currnodes){
			int node_index = iter.first;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
			if (mvpnode){
				ExpandNode(mvpnode, childnodes, node_index);
			}
		}
		ostrm << std::endl;
		currnodes = move(childnodes);
		n+= LPN;
	} while (!currnodes.empty());

	return;
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
size_t mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::MemoryUsage()const{
	std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
	if (m_top != NULL) currnodes[0] = m_top;

	int n_internal=0, n_leaf=0;
	do {
		for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
			int index = iter->first;
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *node = iter->second;

			if (typeid(*node).hash_code() == typeid(MVPInternal<T,BF,PL,LC,LPN,FO,NS>).hash_code())
				n_internal++;
			else
				n_leaf++;
			ExpandNode(node, childnodes, index);
		}
		currnodes = move(childnodes);
	} while (!currnodes.empty());
	
	return  n_internal*sizeof(MVPInternal<T,BF,PL,LC,LPN,FO,NS>)
		+ n_leaf*sizeof(MVPLeaf<T,BF,PL,LC,LPN,FO,NS>)
		+ sizeof(mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,BF,PL,LC,LPN,FO,NS>::SetSync(int n){
	n_sync = n;
}


#endif /* _MVPTREE_H */
