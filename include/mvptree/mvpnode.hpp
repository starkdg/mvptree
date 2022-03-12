#ifndef _MVPNODE_H
#define _MVPNODE_H

#include <map>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <random>
#include "mvptree/mvptree.hpp"
#include "mvptree/datapoint.hpp"

using namespace std;

template<int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
class MVPNode {
public:
	MVPNode(){}

	virtual ~MVPNode(){};

	static MVPNode<BF,PL,LC,LPN,FO,NS>* CreateNode(vector<DataPoint*> &points,
							   map<int,vector<DataPoint*>*> &childpoints,
							   int level, int index);
	
	virtual MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint*> &points,
								   map<int,vector<DataPoint*>*> &childpoints,
								   const int level, const int index) = 0;

	virtual const int GetCount()const = 0;

	virtual void SetChildNode(const int n, MVPNode *node) = 0;

	virtual MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(int n)const = 0;

	virtual const vector<uint64_t> GetVantagePoints()const = 0;
	
	virtual const vector<DataPoint*> GetDataPoints()const = 0;

	virtual const vector<DataPoint*> FilterDataPoints(const uint64_t target_value, const double radius)const = 0;

	virtual void TraverseNode(const uint64_t target_value,
							  const double radius,
							  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
							  const int index,
							  vector<DataPoint*> &results)const = 0;

	virtual const vector<DataPoint*> PurgeDataPoints()=0;

};

template<int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
class MVPInternal : public MVPNode<BF,PL,LC,LPN,FO,NS> {
private:
	int m_nvps;
	uint64_t m_vps[LPN];
	MVPNode<BF,PL,LC,LPN,FO,NS>* m_childnodes[FO];
	double m_splits[LPN][NS];

	void SelectVantagePoints(vector<DataPoint*> &points);
	
	void CalcSplitPoints(const vector<double> &dists, int n, int split_index);

	vector<double> CalcPointDistances(const int vp_index, vector<DataPoint*> &points);

	vector<DataPoint*>* CullPoints(vector<DataPoint*> &list, vector<double> &dists,
								  double split, bool less);
		
	void CollatePoints(vector<DataPoint*> &points,
					   map<int, vector<DataPoint*>*> &childpoints,
					   const int level, const int index);
	
public:
	MVPInternal();
	~MVPInternal(){};
	MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint*> &points,
						   map<int,vector<DataPoint*>*> &childpoints,
						   const int level, const int index);

	const int GetCount()const;

	void SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS> *node);

	MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

	const vector<uint64_t> GetVantagePoints()const;
	
	const vector<DataPoint*> GetDataPoints()const;

	const vector<DataPoint*> FilterDataPoints(const uint64_t target_value, const double radius)const;

	void TraverseNode(const uint64_t target_value, const double radius,
							  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
							  const int index,
							  vector<DataPoint*> &results)const;

	const vector<DataPoint*> PurgeDataPoints();
};

template<int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
class MVPLeaf : public MVPNode<BF,PL,LC,LPN,FO,NS> {
private:
	int m_nvps;
	uint64_t m_vps[PL];
	double m_pdists[PL][LC];

	int m_npoints;
    DataPoint* m_points[LC];

	void SelectVantagePoints(vector<DataPoint*> &points);

	void MarkLeafDistances(vector<DataPoint*> &points);
	
public:
	MVPLeaf();
	~MVPLeaf(){};
	MVPNode<BF,PL,LC,LPN,FO,NS>* AddDataPoints(vector<DataPoint*> &points,
						   map<int,vector<DataPoint*>*> &childpoints,
						   const int level, const int index);

	const int GetCount()const;

	void SetChildNode(const int n, MVPNode<BF,PL,LC,LPN,FO,NS> *node);

	MVPNode<BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

	const vector<uint64_t> GetVantagePoints()const;

	const vector<DataPoint*> GetDataPoints()const;

	const vector<DataPoint*> FilterDataPoints(const uint64_t target_value, const double radius)const;

	void TraverseNode(const uint64_t target_value,const double radius,
					  map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
					  const int index,
					  vector<DataPoint*> &results)const;

	const vector<DataPoint*> PurgeDataPoints();
};


/**
 *  MVPNode implementation 
 *
 **/

bool CompareDistance(const double a, const double b, const bool less){
	if (less) return (a <= b);
	return (a > b);
}

double hamming_distance(const uint64_t a, const uint64_t b){
	return __builtin_popcountll(a^b);
}

/********** MVPNode methods *********************/

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPNode<BF,PL,LC,LPN,FO,NS>::CreateNode(vector<DataPoint*> &points,
							 map<int, vector<DataPoint*>*> &childpoints,
							 int level, int index){
	MVPNode<BF,PL,LC,LPN,FO,NS> *node = NULL;
	if (points.size() <= LC){
		node = new MVPLeaf<BF,PL,LC,LPN,FO,NS>();
	} else {
		node = new MVPInternal<BF,PL,LC,LPN,FO,NS>();
	}

	node = node->AddDataPoints(points, childpoints, level, index);
	if (node == NULL) throw runtime_error("unable to create node");
	return node;
}

/********** MVPInternal methods *******************/
template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPInternal<BF,PL,LC,LPN,FO,NS>::MVPInternal(){
	m_nvps = 0;
	for (int i=0;i < FO;i++) m_childnodes[i] = NULL;
	for (int i=0;i < LPN;i++){
		for (int j=0;j < NS;j++){
			m_splits[i][j] = -1.0;
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(vector<DataPoint*> &points){
	static random_device rd;
	static mt19937 gen(rd());
	static uniform_int_distribution<uint64_t> distrib(0);

	while (m_nvps < LPN){
		uint64_t vp = distrib(gen);
		m_vps[m_nvps++] = vp;
		if (m_nvps < LPN){
			m_vps[m_nvps++] = ~vp;
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
vector<double> MVPInternal<BF,PL,LC,LPN,FO,NS>::CalcPointDistances(const int vp_index, vector<DataPoint*> &points){
	vector<double> results;
	for (DataPoint *dp : points){
		results.push_back(dp->distance(m_vps[vp_index]));
	}
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
vector<DataPoint*>* MVPInternal<BF,PL,LC,LPN,FO,NS>::CullPoints(vector<DataPoint*> &list, vector<double> &dists,
							  double split, bool less){
	vector<DataPoint*> *results = new vector<DataPoint*>();

	vector<DataPoint*>::iterator list_iter = list.begin();
	vector<double>::iterator dist_iter = dists.begin(); 
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
void MVPInternal<BF,PL,LC,LPN,FO,NS>::CollatePoints(vector<DataPoint*> &points,
								map<int, vector<DataPoint*>*> &childpoints,
								const int level, const int index){
	map<int, vector<DataPoint*>*> pnts, pnts2;
	pnts[0] = &points;

	int lengthM = BF - 1;
	int n = 0;
	do {
		for (auto iter=pnts.begin();iter!=pnts.end();iter++){
			int node_index = iter->first;
			vector<DataPoint*> *list = iter->second;

			//int list_size = list->size();

			vector<double> dists = CalcPointDistances(n, *list);
			if (dists.size() > 0){
				CalcSplitPoints(dists, n, node_index);

				double m;
				vector<DataPoint*> *culledpts = NULL;
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
		vector<DataPoint*> *list = iter->second;
		if (list != NULL)
			childpoints[index*FO+i] = list;  
	}
	
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPInternal<BF,PL,LC,LPN,FO,NS>::AddDataPoints(vector<DataPoint*> &points,
									map<int,vector<DataPoint*>*> &childpoints,
									const int level, const int index){
	SelectVantagePoints(points);
	if (m_nvps < LPN) throw invalid_argument("too few points for internal node");
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
const vector<uint64_t> MVPInternal<BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	vector<uint64_t> results;
	for (int i=0;i<m_nvps;i++) results.push_back(m_vps[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPInternal<BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	vector<DataPoint*> results;
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPInternal<BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const uint64_t target_value, const double radius)const{
	vector<DataPoint*> results;
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPInternal<BF,PL,LC,LPN,FO,NS>::TraverseNode(const uint64_t target_value, const double radius,
												   map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
												   const int index, vector<DataPoint*> &results)const{
	int lengthM = BF - 1;
	int n = 0;
	bool *currnodes  = new bool[1];
	currnodes[0] = true;
	do {
		int n_nodes = pow(BF, n);
		int n_childnodes = pow(BF, n+1);
		bool *nextnodes = new bool[n_childnodes];
		for (int i=0;i<n_childnodes;i++) nextnodes[i] = false;
		
		double d =  hamming_distance(target_value, m_vps[n]);

		//int lengthMn = lengthM*n_nodes;
		for (int node_index=0;node_index<n_nodes;node_index++){
			if (currnodes[node_index]){
				if (m_splits[n][node_index*lengthM] >= 0){
					double m = m_splits[n][node_index*lengthM];
					for (int j=0;j<lengthM;j++){
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
			if (child != NULL)
				childnodes[index*FO+i] = child;
		}
	}

	delete[] currnodes;

}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPInternal<BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	vector<DataPoint*> results;
	return results;
}


/********* MVPLeaf methods **********************/

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPLeaf<BF,PL,LC,LPN,FO,NS>::MVPLeaf(){
	m_npoints = m_nvps = 0;
	for (int i=0;i<PL;i++)
		for (int j=0;j<LC;j++)
			m_pdists[i][j] = -1.0;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(vector<DataPoint*> &points){
	static random_device rd;
	static mt19937 gen(rd());
	static uniform_int_distribution<uint64_t> distrib(0);

	while (m_nvps < PL){
		uint64_t vp = distrib(gen);
		m_vps[m_nvps++] = vp;
		if (m_nvps < PL){
			m_vps[m_nvps++] = ~vp;
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::MarkLeafDistances(vector<DataPoint*> &points){
	if (m_npoints + points.size() > LC)
		throw invalid_argument("no. points exceed leaf capacity");
	
	for (int m = 0; m < m_nvps;m++){
		int index = m_npoints;
		for (DataPoint *dp : points){
			m_pdists[m][index++] = dp->distance(m_vps[m]);   
		}
	}
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
MVPNode<BF,PL,LC,LPN,FO,NS>* MVPLeaf<BF,PL,LC,LPN,FO,NS>::AddDataPoints(vector<DataPoint*> &points,
								map<int,vector<DataPoint*>*> &childpoints,
								const int level, const int index){
	SelectVantagePoints(points);
	MVPNode<BF,PL,LC,LPN,FO,NS> *retnode = this;
	if (m_npoints + points.size() <= LC){ 
		// add points to existing leaf
		MarkLeafDistances(points);
		for (DataPoint *dp : points){
			m_points[m_npoints++] = dp;
		}
		points.clear();
	} else {  // create new internal node 

		// get existing points, purge inactive poins
		vector<DataPoint*> pts = PurgeDataPoints();

		// merge points
		for (DataPoint *dp : pts) points.push_back(dp);

		if (points.size() <= LC){
			// clear out points
			m_npoints = 0;
			MarkLeafDistances(points);
			for (DataPoint *dp : points){
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
const vector<uint64_t> MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	vector<uint64_t> results;
	for (int i=0;i<m_nvps;i++) results.push_back(m_vps[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	vector<DataPoint*> results;
	for (int i=0;i<m_npoints;i++) results.push_back(m_points[i]);
	return results;
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const uint64_t target_value, const double radius)const{
	vector<DataPoint*> results;

	double qdists[PL];
	for (int i=0;i<m_nvps;i++){
		qdists[i] = hamming_distance(m_vps[i], target_value);
	}
	
	for (int j=0;j < (int)m_npoints;j++){
		bool skip = false;
		if (!(m_points[j]->IsActive())) continue;
		
		for (int i=0;i<m_nvps;i++){
			if (!(m_pdists[i][j] >= qdists[i] - radius) && (m_pdists[i][j] <= qdists[i] + radius)){
				skip = true;
				break;
			}
		}
		if (!skip){
			if (m_points[j]->distance(target_value) <= radius){
				results.push_back(m_points[j]);
			}
		}
	}
	
	return results;
}


template<int BF,int PL,int LC,int LPN,int FO,int NS>
void MVPLeaf<BF,PL,LC,LPN,FO,NS>::TraverseNode(const uint64_t target_value, const double radius,
						   map<int, MVPNode<BF,PL,LC,LPN,FO,NS>*> &childnodes,
						   const int index, vector<DataPoint*> &results)const{
	vector<DataPoint*> fnd = FilterDataPoints(target_value, radius);
	for (DataPoint* dp : fnd) results.push_back(dp);
}

template<int BF,int PL,int LC,int LPN,int FO,int NS>
const vector<DataPoint*> MVPLeaf<BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	vector<DataPoint*> results;
	for (int i=0;i<m_npoints;i++){
		if (!(m_points[i]->IsActive()))
			delete m_points[i];
		else
			results.push_back(m_points[i]);
	}
	return results;
}


#endif /* _MVPNODE_H */
