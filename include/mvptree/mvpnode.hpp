#ifndef _MVPNODE_H
#define _MVPNODE_H

#include <map>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "mvptree/mvptree.hpp"
#include "mvptree/datapoint.hpp"

using namespace std;

template<int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
class MVPNode {
protected:

	int m_nvps;
	DataPoint<PL>* m_vps[LPN];

	void SelectVantagePoints(vector<DataPoint<PL>*> &points);
	
	void MarkPointDistances(vector<DataPoint<PL>*> &points, const int level);

public:
	MVPNode():m_nvps(0){}

	virtual ~MVPNode(){};

	static MVPNode<BF,PL,LC,LPN,FO,NS>* CreateNode(vector<DataPoint<PL>*> &points,
							   map<int,vector<DataPoint<PL>*>*> &childpoints,
							   int level, int index);
	
	virtual MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint<PL>*> &points,
								   map<int,vector<DataPoint<PL>*>*> &childpoints,
								   const int level, const int index) = 0;

	virtual const int GetCount()const = 0;

	virtual void SetChildNode(const int n, MVPNode *node) = 0;

	virtual MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(int n)const = 0;

	virtual const vector<DataPoint<PL>*> GetVantagePoints()const = 0;
	
	virtual const vector<DataPoint<PL>*> GetDataPoints()const = 0;

	virtual void FilterDataPoints(const DataPoint<PL> &target,
										vector<double> &tpath,
										const double radius,
										const int level,
										vector<DataPoint<PL>*> &results)const = 0;
 
	virtual void TraverseNode(const DataPoint<PL> &target,
							  const double radius,
							  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
							  vector<double> *tpath,
							  map<int, vector<double>*> &child_tdistances,
							  const int index,
							  const int level,
							  vector<DataPoint<PL>*> &results)const = 0;

	virtual const vector<DataPoint<PL>*> PurgeDataPoints()=0;

};

template<int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
class MVPInternal : public MVPNode<BF,PL,LC,LPN,FO,NS> {
private:

	MVPNode<BF,PL,LC,LPN,FO,NS>* m_childnodes[FO];
	double m_splits[LPN][NS];

	void CalcSplitPoints(const vector<double> &dists, int n, int split_index);

	vector<double> CalcPointDistances(DataPoint<PL> &vp, vector<DataPoint<PL>*> &points);

	vector<DataPoint<PL>*>* CullPoints(vector<DataPoint<PL>*> &list, vector<double> &dists,
								  double split, bool less);
		
	void CollatePoints(vector<DataPoint<PL>*> &points,
					   map<int, vector<DataPoint<PL>*>*> &childpoints,
					   const int level, const int index);
	
public:
	MVPInternal();
	~MVPInternal(){};
	MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint<PL>*> &points,
						   map<int,vector<DataPoint<PL>*>*> &childpoints,
						   const int level, const int index);

	const int GetCount()const;

	void SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS> *node);

	MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

	const vector<DataPoint<PL>*> GetVantagePoints()const;
	
	const vector<DataPoint<PL>*> GetDataPoints()const;

	void FilterDataPoints(const DataPoint<PL> &target,
												  vector<double> &tpath,
												  const double radius,
												  const int level,
												  vector<DataPoint<PL>*> &results)const;

	void TraverseNode(const DataPoint<PL> &target,const double radius,
					  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
					  vector<double> *tpath,
					  map<int, vector<double>*> &child_tdistances,
					  const int index,
					  const int level,
					  vector<DataPoint<PL>*> &results)const;

	const vector<DataPoint<PL>*> PurgeDataPoints();
};

template<int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
class MVPLeaf : public MVPNode<BF,PL,LC,LPN,FO,NS> {
private:

	double m_pdists[LPN][LC];
	int m_npoints;
    DataPoint<PL>* m_points[LC];

	void MarkLeafDistances(vector<DataPoint<PL>*> &points);
	
public:
	MVPLeaf();
	~MVPLeaf(){};
	MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint<PL>*> &points,
						   map<int,vector<DataPoint<PL>*>*> &childpoints,
						   const int level, const int index);

	const int GetCount()const;

	void SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS> *node);

	MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

	const vector<DataPoint<PL>*> GetVantagePoints()const;

	const vector<DataPoint<PL>*> GetDataPoints()const;

	void FilterDataPoints(const DataPoint<PL> &target,
												  vector<double> &tpath,
												  const double radius,
												  const int level,
												  vector<DataPoint<PL>*> &results)const;

	void TraverseNode(const DataPoint<PL> &target,const double radius,
					  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
					  vector<double> *tpath,
					  map<int, vector<double>*> &child_tdistances,
 					  const int index,
					  const int level,
					  vector<DataPoint<PL>*> &results)const;

	const vector<DataPoint<PL>*> PurgeDataPoints();
};


/**
 *  MVPNode implementation 
 *
 **/

bool CompareDistance(const double a, const double b, const bool less){
	if (less) return (a <= b);
	return (a > b);
}

/********** MVPNode methods *********************/

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPNode<BF,PL,LC,LPN,FO,NS>::CreateNode(vector<DataPoint<PL>*> &points,
							 map<int, vector<DataPoint<PL>*>*> &childpoints,
							 int level, int index){
	MVPNode<BF,PL,LC,LPN,FO,NS> *node = NULL;
	if (points.size() <= LC + LPN){
		node = new MVPLeaf<BF,PL,LC,LPN,FO,NS>();
	} else {
		node = new MVPInternal<BF,PL,LC,LPN,FO,NS>();
	}

	node = node->AddDataPoints(points, childpoints, level, index);
	if (node == NULL) throw runtime_error("unable to create node");
	return node;
}


template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPNode<BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(vector<DataPoint<PL>*> &points){
	while (this->m_nvps < LPN && points.size() > 0){
		this->m_vps[this->m_nvps++] = points.back();
		points.pop_back();
	}
}


/**
// alternate select vantage points
template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPNode<BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(vector<DataPoint<PL>*> &points){
	const int limit = 50;
	if (this->m_nvps < LPN && points.size() > 0){

		DataPoint<PL> *vp = points.back();
		
		this->m_vps[this->m_nvps++] = vp;
		points.pop_back();

		while (this->m_nvps < LPN && points.size() > 0) {
			double max_distance = -1.0;
			int max_pos = 0;
			int nlimit = ((int)points.size() <= limit) ? (int)points.size() : limit;
			for (int i=0;i < nlimit;i++){
				double d = vp->distance(points[i]);
				if (d > max_distance){
					max_distance = d;
					max_pos = i;
				}
			}

			vp = points[max_pos];
			points[max_pos] = points.back();
			points.pop_back();
			this->m_vps[m_nvps++] = vp;
		}
	}
}
**/

template<int BF, int PL, int LC, int LPN, int FO, int NS>
void MVPNode<BF,PL,LC,LPN,FO,NS>::MarkPointDistances(vector<DataPoint<PL>*> &points, const int level){
	for (int i=0;i < this->m_nvps;i++){
		if (level + i < PL){
			for (DataPoint<PL> *pnt : points){
				double d = pnt->distance(m_vps[i]);
				pnt->SetPath(d, level+i);
			}
		}
	}
}

/********** MVPInternal methods *******************/
template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPInternal<BF,PL,LC,LPN,FO,NS>::MVPInternal(){
	for (int i=0;i < FO;i++) m_childnodes[i] = NULL;
	for (int i=0;i < LPN;i++){
		for (int j=0;j < NS;j++){
			m_splits[i][j] = -1.0;
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::CalcSplitPoints(const vector<double> &dists, int n, int split_index){
	int lengthM = BF - 1;
	if (dists.size() > 0){
		if (m_splits[n][split_index*lengthM] == -1){
			vector<double> tmpdists = dists;
			sort(tmpdists.begin(), tmpdists.end());
			double factor = (double)tmpdists.size()/(double)BF;
			for (int i=0;i<lengthM;i++){
				double pos = (i+1)*factor;
				int lo = floor(pos);
				int hi = (pos <= tmpdists.size()-1) ? ceil(pos) : 0;
				m_splits[n][split_index*lengthM+i] = (tmpdists[lo] + tmpdists[hi])/2.0;
			}
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
vector<double> MVPInternal<BF,PL,LC,LPN,FO,NS>::CalcPointDistances(DataPoint<PL> &vp, vector<DataPoint<PL>*> &points){
	vector<double> results;
	for (DataPoint<PL> *dp : points){
		results.push_back(vp.distance(dp));
	}
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
vector<DataPoint<PL>*>* MVPInternal<BF,PL,LC,LPN,FO,NS>::CullPoints(vector<DataPoint<PL>*> &list, vector<double> &dists,
							  double split, bool less){
	vector<DataPoint<PL>*> *results = new vector<DataPoint<PL>*>();

	auto list_iter = list.begin();
	auto dist_iter = dists.begin(); 
	while (list_iter != list.end() && dist_iter != dists.end()){
		if (CompareDistance(*dist_iter,split,less)){
			results->push_back(*list_iter);
			list_iter = list.erase(list_iter);
			dist_iter = dists.erase(dist_iter);
		} else {
			list_iter++;
			dist_iter++;
		}
	}
	if (results->size() > 0){
		return results;
	}
	delete results;
	return NULL;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::CollatePoints(vector<DataPoint<PL>*> &points,
								map<int, vector<DataPoint<PL>*>*> &childpoints,
								const int level, const int index){
	map<int, vector<DataPoint<PL>*>*> pnts, pnts2;
	pnts[0] = &points;

	int lengthM = BF - 1;
	int n = 0;
	do {
		for (auto iter=pnts.begin();iter!=pnts.end();iter++){
			int node_index = iter->first;
			vector<DataPoint<PL>*> *list = iter->second;

			//int list_size = list->size();

			vector<double> dists = CalcPointDistances(*(this->m_vps[n]), *list);
			if (dists.size() > 0){
				CalcSplitPoints(dists, n, node_index);

				double m;
				vector<DataPoint<PL>*> *culledpts = NULL;
				for (int j=0;j<lengthM;j++){
					m = m_splits[n][node_index*lengthM+j];

					culledpts = CullPoints(*list, dists, m, true);
					if (culledpts != NULL){
						pnts2[node_index*BF+j] = culledpts;
					}
				}
				m = m_splits[n][node_index*lengthM+lengthM-1];
				culledpts = CullPoints(*list, dists, m, false);
				if (culledpts != NULL){
					pnts2[node_index*BF+BF-1] = culledpts;
				}
			}
			if (list->size() > 0)
				throw length_error("not fully collated");

			if (n > 0) delete list;
		}
		pnts = move(pnts2);
		n++;
	} while (n < LPN);

	for (auto iter=pnts.begin();iter!=pnts.end();iter++){
		int i=iter->first;
		vector<DataPoint<PL>*> *list = iter->second;
		if (list != NULL)
			childpoints[index*FO+i] = list;  
	}
	
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPInternal<BF,PL,LC,LPN,FO,NS>::AddDataPoints(vector<DataPoint<PL>*> &points,
									map<int,vector<DataPoint<PL>*>*> &childpoints,
									const int level, const int index){
	this->SelectVantagePoints(points);
	this->MarkPointDistances(points, level);
	if (this->m_nvps < LPN) throw invalid_argument("too few points for internal node");
	CollatePoints(points, childpoints, level, index);
	points.clear();
	return this;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const int MVPInternal<BF,PL,LC,LPN,FO,NS>::GetCount()const{
	return 0;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS> *node){
	if (n < 0 || n >= FO) throw invalid_argument("index out of range");
	m_childnodes[n] = node;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPInternal<BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{
	if (n < 0 || n >= FO) throw invalid_argument("index out of range");
	return m_childnodes[n];
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPInternal<BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	vector<DataPoint<PL>*> results;
	for (int i=0;i < this->m_nvps;i++) results.push_back(this->m_vps[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPInternal<BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	vector<DataPoint<PL>*> results;
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const DataPoint<PL> &target,
													   vector<double> &tpath,
													   const double radius,
													   const int level,
													   vector<DataPoint<PL>*> &results)const{
	for (int i=0;i < this->m_nvps;i++){
		double d = target.distance(this->m_vps[i]);
		tpath.push_back(d);
		if (this->m_vps[i]->IsActive() && d <= radius)
			results.push_back(this->m_vps[i]);
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::TraverseNode(const DataPoint<PL> &target, const double radius,
												   map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
												   vector<double> *tpath,
												   map<int, vector<double>*> &child_tdistances,
												   const int index,
												   const int level,
												   vector<DataPoint<PL>*> &results)const{
	int lengthM = BF - 1;
	int n = 0;
	bool *currnodes  = new bool[1];
	currnodes[0] = true;


	vector<double> *path = (tpath != NULL) ? tpath : new vector<double>();

	FilterDataPoints(target, *path, radius, level, results);
	
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+1);
		bool *nextnodes = new bool[n_childnodes];

		for (int i=0;i<n_childnodes;i++)
			nextnodes[i] = false;

		double d = path->at(level+n);

		//int lengthMn = lengthM*n_nodes;
		for (int node_index=0;node_index < n_nodes;node_index++){
			if (currnodes[node_index]){
				if (m_splits[n][node_index*lengthM] >= 0){
					double m = m_splits[n][node_index*lengthM];
					for (int j=0;j < lengthM;j++){
						m = m_splits[n][node_index*lengthM+j];
						if (d <= m + radius) nextnodes[node_index*BF+j] = true;
					}
					if (d > m - radius) nextnodes[node_index*BF+BF-1] = true;
				}
			}
		}

		delete[] currnodes;
		currnodes = nextnodes;
		n++;
	} while (n < LPN);

	for (int i=0;i < FO;i++){
		if (currnodes[i]){
			MVPNode<BF,PL,LC,LPN,FO,NS> *child = GetChildNode(i);
			if (child != NULL) {
				childnodes[index*FO+i] = child;

				vector<double> *tmppath = new vector<double>();
				for (int j=0;j < (int)path->size();j++){
					tmppath->push_back(path->at(j));
				}

				child_tdistances[index*FO+i] = tmppath;
			}
		}
	   
	}

	delete[] currnodes;
	delete path;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPInternal<BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	vector<DataPoint<PL>*> results;
	for (int i=0;i < this->m_nvps;i++){
		if (!(this->m_vps[i]->IsActive()))
			delete this->m_vps[i];
		else
			results.push_back(this->m_vps[i]);
	}
	
	return results;
}


/********* MVPLeaf methods **********************/

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPLeaf<BF,PL,LC,LPN,FO,NS>::MVPLeaf():m_npoints(0){
	for (int i=0;i < LPN;i++)
		for (int j=0;j<LC;j++)
			m_pdists[i][j] = -1.0;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::MarkLeafDistances(vector<DataPoint<PL>*> &points){
	if (m_npoints + points.size() > LC)
		throw invalid_argument("no. points exceed leaf capacity");
	
	for (int m = 0; m < this->m_nvps;m++){
		int index = m_npoints;
		for (DataPoint<PL> *dp : points){
			m_pdists[m][index++] = this->m_vps[m]->distance(dp);  
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPLeaf<BF,PL,LC,LPN,FO,NS>::AddDataPoints(vector<DataPoint<PL>*> &points,
								map<int,vector<DataPoint<PL>*>*> &childpoints,
								const int level, const int index){
	this->SelectVantagePoints(points);
	MVPNode<BF,PL,LC,LPN,FO,NS> *retnode = this;
	if (m_npoints + points.size() <= LC){ 
		// add points to existing leaf
		MarkLeafDistances(points);
		this->MarkPointDistances(points, level);
		for (DataPoint<PL> *dp : points){
			m_points[m_npoints++] = dp;
		}
		points.clear();
	} else {  // create new internal node 

		// get existing points, purge inactive poins
		vector<DataPoint<PL>*> pts = PurgeDataPoints();

		// merge points
		for (DataPoint<PL> *dp : pts) points.push_back(dp);

		if (points.size() <= LPN + LC){
			// clear out points
			this->m_npoints = 0;
			this->m_nvps = 0;
			this->SelectVantagePoints(points);
			MarkLeafDistances(points);
			this->MarkPointDistances(points, level);
			for (DataPoint<PL> *dp : points){
				m_points[m_npoints++] = dp;
			}
			points.clear();
		} else {
			retnode = new MVPInternal<BF,PL,LC,LPN,FO,NS>();
			retnode = retnode->AddDataPoints(points, childpoints, level, index);
		}
	}
	if (retnode == NULL) throw runtime_error("unable to create node");
	
	return retnode;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const int MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetCount()const{
	return m_npoints;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS>* node){}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{return NULL;}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	vector<DataPoint<PL>*> results;
	for (int i=0;i < this->m_nvps;i++) results.push_back(this->m_vps[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	vector<DataPoint<PL>*> results;
	for (int i=0;i<m_npoints;i++) results.push_back(m_points[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const DataPoint<PL> &target,
												   vector<double> &tpath,
												   const double radius,
												   const int level,
												   vector<DataPoint<PL>*> &results)const{
	for (int i=0;i < this->m_nvps;i++){
		double d = target.distance(this->m_vps[i]);
		tpath.push_back(d);
		if (this->m_vps[i]->IsActive() && d <= radius)
			results.push_back(this->m_vps[i]);
	}

	int pathlimit = (tpath.size() <= PL) ? tpath.size() : PL;
   	for (int index=0;index < (int)m_npoints;index++){
		bool skip = false;
		DataPoint<PL> *dp = m_points[index];

		if (dp->IsActive()){
			// filter using precomputed distances in node
			for (int i=0;i < this->m_nvps;i++){
				if (!(m_pdists[i][index] >=  tpath[level+i] - radius) && (m_pdists[i][index] <= tpath[level+i] + radius)){
					skip = true;
					break;
				}
			}

			if (!skip){
				// filter using precomputed path distances between target and vantage points
				// from top down to current place in tree
				for (int i=0;i < pathlimit;i++){
					if (!(dp->GetPath(i) >= tpath[i] - radius && dp->GetPath(i) <= tpath[i] + radius)){
						skip = true;
						break;
					}
				}
			}

			if (!skip){
				// still not ruled out
				if (target.distance(dp) <= radius){
					results.push_back(dp);
				}
			}
		}
	}
}


template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::TraverseNode(const DataPoint<PL> &target, const double radius,
											   map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
											   vector<double> *tpath,
											   map<int, vector<double>*> &child_tdistances,
											   const int index,
											   const int level,
											   vector<DataPoint<PL>*> &results)const{

	vector<double> *path = (tpath) ? tpath : new vector<double>();
	FilterDataPoints(target, *path, radius, level, results);
	delete path;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint<PL>*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	vector<DataPoint<PL>*> results;
	for (int i=0;i < this->m_nvps;i++){
		if (!(this->m_vps[i]->IsActive()))
			delete this->m_vps[i];
		else
			results.push_back(this->m_vps[i]);
	}
	for (int i=0;i<m_npoints;i++){
		if (!(m_points[i]->IsActive()))
			delete m_points[i];
		else
			results.push_back(m_points[i]);
	}
	return results;
}



#endif /* _MVPNODE_H */
