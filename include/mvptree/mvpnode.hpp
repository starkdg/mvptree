/**
    MVPTree
    Copyright (C) 2022  David G. Starkweather

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
**/
#ifndef _MVPNODE_H
#define _MVPNODE_H

#include <map>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "mvptree/datapoint.hpp"


namespace mvp {

	template<typename T, int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
	class MVPNode {
	protected:

		int m_nvps;
		std::array<vp_t<T>,LPN>  m_vps;

		void SelectVantagePoints(std::vector<datapoint_t<T,PL>> &points);
		
		void MarkPointDistances(std::vector<datapoint_t<T,PL>> &points, const int level);

	public:
		MVPNode():m_nvps(0){}

		virtual ~MVPNode(){};

		static MVPNode<T,BF,PL,LC,LPN,FO,NS>* CreateNode(std::vector<datapoint_t<T,PL>> &points,
														 std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
														 int level,
														 int index);
	
		virtual MVPNode<T,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
															 std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
															 const int level,
															 const int index) = 0;

		virtual const size_t nbytes()const = 0;

		virtual const int Size()const = 0;

		virtual const bool IsLeaf()const = 0;
	
		virtual void SetChildNode(const int n, MVPNode<T,BF,PL,LC,LPN,FO,NS> *node) = 0;
		
		virtual MVPNode<T,BF,PL,LC,LPN,FO,NS>* GetChildNode(int n)const = 0;

		virtual const std::vector<vp_t<T>> GetVantagePoints()const = 0;
		
		virtual const std::vector<datapoint_t<T,PL>> GetDataPoints()const = 0;

		virtual void FilterDataPoints(const T &target,
									  std::vector<double> &tpath,
									  const double radius,
									  const int level,
									  const bool delete_points,
									  std::vector<item_t<T>> &results) = 0;
 
		virtual void TraverseNode(const T &target,
								  const double radius,
								  std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
								  std::vector<double> *tpath,
								  std::map<int, std::vector<double>*> &child_tdistances,
								  const int index,
								  const int level,
								  const bool delete_points,
								  std::vector<item_t<T>> &results) = 0;

		virtual const std::vector<datapoint_t<T,PL>> PurgeDataPoints()=0;

	};

	template<typename T, int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
	class MVPInternal : public MVPNode<T,BF,PL,LC,LPN,FO,NS> {
	private:

		std::array<MVPNode<T,BF,PL,LC,LPN,FO,NS>*,FO> m_childnodes;

		// LPN x NS;
		std::array<std::array<double,NS>,LPN> m_splits;
	
		void CalcSplitPoints(const std::vector<double> &dists, int n, int split_index);

		std::vector<double> CalcPointDistances(vp_t<T> &vp, std::vector<datapoint_t<T,PL>> &points);

		std::vector<datapoint_t<T,PL>>* CullPoints(std::vector<datapoint_t<T,PL>> &list, std::vector<double> &dists,
												   double split, bool less);
		
		void CollatePoints(std::vector<datapoint_t<T,PL>> &points,
						   std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints,
						   const int level, const int index);
	
	public:
		MVPInternal();
		~MVPInternal(){};
		
		MVPNode<T,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
													 std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
													 const int level,
													 const int index);

		const size_t nbytes()const;
		
		const int Size()const;

		const bool IsLeaf() const;
	
		void SetChildNode(const int n, MVPNode<T,BF,PL,LC,LPN,FO,NS> *node);

		MVPNode<T,BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

		const std::vector<vp_t<T>> GetVantagePoints()const;
	
		const std::vector<datapoint_t<T,PL>> GetDataPoints()const;

		void FilterDataPoints(const T &target,
							  std::vector<double> &tpath,
							  const double radius,
							  const int level,
							  const bool delete_points,
							  std::vector<item_t<T>> &results);

		void TraverseNode(const T &target,const double radius,
						  std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
						  std::vector<double> *tpath,
						  std::map<int, std::vector<double>*> &child_tdistances,
						  const int index,
						  const int level,
						  const bool delete_points,
						  std::vector<item_t<T>> &results);

		const std::vector<datapoint_t<T,PL>> PurgeDataPoints();
	};

	template<typename T,int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
	class MVPLeaf : public MVPNode<T,BF,PL,LC,LPN,FO,NS> {
	private:

		// PN x LC;
		std::array<std::array<double,LC>,LPN> m_pdists;
		int m_npoints;
		std::array<datapoint_t<T,PL>,LC>  m_points;

		void MarkLeafDistances(std::vector<datapoint_t<T,PL>> &points);
	
	public:
		MVPLeaf();
		~MVPLeaf(){};
		MVPNode<T,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
													 std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
													 const int level, const int index);

		const size_t nbytes()const;

		const int Size()const;

		const bool IsLeaf()const;
	
		void SetChildNode(const int n, MVPNode<T,BF,PL,LC,LPN,FO,NS> *node);
		
		MVPNode<T,BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

		const std::vector<vp_t<T>> GetVantagePoints()const;
		
		const std::vector<datapoint_t<T,PL>> GetDataPoints()const;

		void FilterDataPoints(const T &target,
							  std::vector<double> &tpath,
							  const double radius,
							  const int level,
							  const bool delete_points,
							  std::vector<item_t<T>> &results);

		void TraverseNode(const T &target,const double radius,
						  std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
						  std::vector<double> *tpath,
						  std::map<int, std::vector<double>*> &child_tdistances,
						  const int index,
						  const int level,
						  const bool delete_points,
						  std::vector<item_t<T>> &results);

		const std::vector<datapoint_t<T,PL>> PurgeDataPoints();
	};


	/**
	 *  MVPNode implementation 
	 *
	 **/

	bool CompareDistance(const double a, const double b, const bool less){
		return (less) ? (a <= b) : (a > b);
	}

}
	
/********** MVPNode methods *********************/

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>::CreateNode(std::vector<mvp::datapoint_t<T,PL>> &points,
																				   std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
																				   int level,
																				   int index){
	MVPNode<T,BF,PL,LC,LPN,FO,NS> *node = NULL;
	if (points.size() <= LC + LPN){
		node = new MVPLeaf<T,BF,PL,LC,LPN,FO,NS>();
	} else {
		node = new MVPInternal<T,BF,PL,LC,LPN,FO,NS>();
	}

	node = node->AddDataPoints(points, childpoints, level, index);
	if (node == NULL) throw std::runtime_error("unable to create node");
	return node;
}


template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(std::vector<mvp::datapoint_t<T,PL>> &points){
   
	if (this->m_nvps == 0 && points.size() > 0){
		this->m_vps[this->m_nvps++] = { points.back().id, points.back().key };
		points.pop_back();
	}

	int slimit = 10;
	while (this->m_nvps < LPN && points.size() > 0){
		int max_pos = 0;
		double maxd = 0;
		slimit = ((int)points.size() <= slimit) ? points.size() : slimit;
		for (int i=0;i < slimit;i++){
			double d = points[i].distance(this->m_vps[this->m_nvps-1]);
			if (d > maxd){
				max_pos = i;
				maxd = d;
			}
		}
		this->m_vps[this->m_nvps++] = { points[max_pos].id, points[max_pos].key };
		points[points.size()-1] = points[max_pos];
		points.pop_back();
	}
}

template<typename T,int BF, int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>::MarkPointDistances(std::vector<mvp::datapoint_t<T,PL>> &points, const int level){
	for (int i=0;i < this->m_nvps;i++){
		if (level + i < PL){
			for (mvp::datapoint_t<T,PL> &pnt : points){
				pnt.dists[level+i] = pnt.distance(m_vps[i]); 
			}
		}
	}
}

/********** MVPInternal methods *******************/
template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::MVPInternal(){
	for (int i=0;i < FO;i++) m_childnodes[i] = NULL;
	for (int i=0;i < LPN;i++){
		for (int j=0;j < NS;j++){
			m_splits[i][j] = -1.0;
		}
	}
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::CalcSplitPoints(const std::vector<double> &dists, int n, int split_index){
	int lengthM = BF - 1;
	if (dists.size() > 0){
		if (m_splits[n][split_index*lengthM] == -1){
			std::vector<double> tmpdists = dists;
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

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
std::vector<double> mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::CalcPointDistances(vp_t<T> &vp, std::vector<mvp::datapoint_t<T,PL>> &points){
	std::vector<double> results;
	for (mvp::datapoint_t<T,PL> &dp : points){
		results.push_back(dp.distance(vp));
	}
	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
std::vector<mvp::datapoint_t<T,PL>>* mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::CullPoints(std::vector<mvp::datapoint_t<T,PL>> &list,
																		 std::vector<double> &dists,
																		 double split,
																		 bool less){
	std::vector<mvp::datapoint_t<T,PL>> *results = new std::vector<mvp::datapoint_t<T,PL>>();

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

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::CollatePoints(std::vector<mvp::datapoint_t<T,PL>> &points,
								std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
								const int level, const int index){
	std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> pnts, pnts2;
	pnts[0] = &points;

	int lengthM = BF - 1;
	int n = 0;
	do {
		for (auto iter=pnts.begin();iter!=pnts.end();iter++){
			int node_index = iter->first;
			std::vector<mvp::datapoint_t<T,PL>> *list = iter->second;

			//int list_size = list->size();

			std::vector<double> dists = CalcPointDistances((this->m_vps[n]), *list);
			if (dists.size() > 0){
				CalcSplitPoints(dists, n, node_index);

				double m;
				std::vector<mvp::datapoint_t<T,PL>> *culledpts = NULL;
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
				throw std::length_error("not fully collated");

			if (n > 0) delete list;
		}
		pnts = move(pnts2);
		n++;
	} while (n < LPN);

	for (auto iter=pnts.begin();iter!=pnts.end();iter++){
		int i=iter->first;
		std::vector<mvp::datapoint_t<T,PL>> *list = iter->second;
		if (list != NULL)
			childpoints[index*FO+i] = list;  
	}
	
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::AddDataPoints(std::vector<mvp::datapoint_t<T,PL>> &points,
																				std::map<int,std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
																				const int level,
																				const int index){
	this->SelectVantagePoints(points);
	this->MarkPointDistances(points, level);
	if (this->m_nvps < LPN) throw std::invalid_argument("too few points for internal node");
	CollatePoints(points, childpoints, level, index);
	points.clear();
	return this;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const size_t mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::nbytes()const{
	return sizeof(MVPInternal<T,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::Size()const{
	int total = 0;
	for (int i=0;i < this->m_nvps;i++){
		if (this->m_vps[i].active) total++;
	}
	return total;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const bool mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::IsLeaf()const{
	return false;
}


template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<T,BF,PL,LC,LPN,FO,NS> *node){
	if (n < 0 || n >= FO) throw std::invalid_argument("index out of range");
	m_childnodes[n] = node;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{
	if (n < 0 || n >= FO) throw std::invalid_argument("index out of range");
	return m_childnodes[n];
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::vp_t<T>> mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	std::vector<mvp::vp_t<T>> results;
	for (int i=0;i < this->m_nvps;i++) results.push_back(this->m_vps[i]);
	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	std::vector<mvp::datapoint_t<T,PL>> results;
	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const T &target,
														 std::vector<double> &tpath,
														 const double radius,
														 const int level,
														 const bool delete_points,
														 std::vector<item_t<T>> &results){
	for (int i=0;i < this->m_nvps;i++){
		double d = this->m_vps[i].distance(target); 
		tpath.push_back(d);
		if (this->m_vps[i].active && d <= radius){
			mvp::vp_t<T> vp = this->m_vps[i];
			if (delete_points) this->m_vps[i].active = false;
			results.push_back({vp.id, vp.key});
		}
	}
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::TraverseNode(const T &target,
													 const double radius,
													 std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
													 std::vector<double> *tpath,
													 std::map<int, std::vector<double>*> &child_tdistances,
													 const int index,
													 const int level,
													 const bool delete_points,
													 std::vector<item_t<T>> &results){
	int lengthM = BF - 1;
	int n = 0;
	bool *currnodes  = new bool[1];
	currnodes[0] = true;


	std::vector<double> *path = (tpath != NULL) ? tpath : new std::vector<double>();

	FilterDataPoints(target, *path, radius, level, delete_points, results);
	
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
			MVPNode<T,BF,PL,LC,LPN,FO,NS> *child = GetChildNode(i);
			if (child != NULL) {
				childnodes[index*FO+i] = child;

				std::vector<double> *tmppath = new std::vector<double>();
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

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	std::vector<mvp::datapoint_t<T,PL>> results;
	for (int i=0;i < this->m_nvps;i++){
		if (this->m_vps[i].active){
			vp_t<T> vp = this->m_vps[i];
			results.push_back({ vp.id, vp.key} );
		}
	}

	return results;
}

/********* MVPLeaf methods **********************/

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::MVPLeaf(){
	m_npoints = 0;
	for (int i=0;i < LPN;i++)
		for (int j=0;j<LC;j++)
			m_pdists[i][j] = -1.0;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::MarkLeafDistances(std::vector<mvp::datapoint_t<T,PL>> &points){
	if (m_npoints + points.size() > LC)
		throw std::invalid_argument("no. points exceed leaf capacity");
	
	for (int m = 0; m < this->m_nvps;m++){
		int index = m_npoints;
		vp_t<T> vp = this->m_vps[m];
		for (mvp::datapoint_t<T,PL> &dp : points){
			m_pdists[m][index++] = dp.distance(vp);  
		}
	}
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::AddDataPoints(std::vector<mvp::datapoint_t<T,PL>> &points,
								std::map<int,std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
								const int level, const int index){
	this->SelectVantagePoints(points);
	MVPNode<T,BF,PL,LC,LPN,FO,NS> *retnode = this;
	if (m_npoints + points.size() <= LC){ 
		// add points to existing leaf
		MarkLeafDistances(points);
		this->MarkPointDistances(points, level);
		for (mvp::datapoint_t<T,PL> &dp : points){
			m_points[m_npoints++] = dp;
		}
		points.clear();
	} else {  // create new internal node 

		// get existing points, purge inactive poins
		std::vector<mvp::datapoint_t<T,PL>> pts = PurgeDataPoints();

		// merge points
		for (mvp::datapoint_t<T,PL> &dp : pts)
			points.push_back(dp);

		if (points.size() <= LPN + LC){
			// clear out points
			this->m_npoints = 0;
			this->m_nvps = 0;
			this->SelectVantagePoints(points);
			MarkLeafDistances(points);
			this->MarkPointDistances(points, level);
			for (mvp::datapoint_t<T,PL> &dp : points){
				m_points[m_npoints++] = dp;
			}
			points.clear();
		} else {
			retnode = new mvp::MVPInternal<T,BF,PL,LC,LPN,FO,NS>();
			retnode = retnode->AddDataPoints(points, childpoints, level, index);
		}
	}
	if (retnode == NULL) throw std::runtime_error("unable to create node");
	
	return retnode;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const size_t mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::nbytes()const{
	return sizeof(mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::Size()const{
	int total = 0;
	for (int i=0;i < this->m_nvps;i++){
		if (this->m_vps[i].active) total++;
	}
	for (int i=0;i < this->m_npoints;i++){
		if (m_points[i].active) total++;
	}

	return total;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const bool mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::IsLeaf()const{
	return true;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<T,BF,PL,LC,LPN,FO,NS>* node){}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,BF,PL,LC,LPN,FO,NS>* mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{return NULL;}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::vp_t<T>> mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
	std::vector<vp_t<T>> results;
	for (int i=0;i < this->m_nvps;i++)
		results.push_back(this->m_vps[i]);

	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
	std::vector<mvp::datapoint_t<T,PL>> results;
	for (int i=0;i<m_npoints;i++)
		results.push_back(m_points[i]);
	return results;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const T &target,
													 std::vector<double> &tpath,
													 const double radius,
													 const int level,
													 const bool delete_points,
													 std::vector<item_t<T>> &results){
	for (int i=0;i < this->m_nvps;i++){
		double d = this->m_vps[i].distance(target);
		tpath.push_back(d);
		if (this->m_vps[i].active && d <= radius){
			vp_t<T> vp = this->m_vps[i];
			if (delete_points){
				vp.active = false;
			}
			results.push_back({vp.id, vp.key });
		}
	}

	int pathlimit = (tpath.size() <= PL) ? tpath.size() : PL;
   	for (int index=0;index < (int)m_npoints;index++){
		bool skip = false;
		mvp::datapoint_t<T,PL> dp = m_points[index];

		if (dp.active){
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
					if (!(dp.dists[i] >= tpath[i] - radius && dp.dists[i] <= tpath[i] + radius)){
						skip = true;
						break;
					}
				}
			}

			if (!skip){
				// still not ruled out
				if (dp.distance(target) <= radius){
					if (delete_points) {
						m_points[index].active = false;
					}
					results.push_back({ dp.id, dp.key });
				}
			}
		}
	}
}


template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::TraverseNode(const T &target,
												 const double radius,
												 std::map<int, MVPNode<T,BF,PL,LC,LPN,FO,NS>*> &childnodes,
												 std::vector<double> *tpath,
												 std::map<int, std::vector<double>*> &child_tdistances,
												 const int index,
												 const int level,
												 const bool delete_points,
												 std::vector<item_t<T>> &results){

	std::vector<double> *path = (tpath) ? tpath : new std::vector<double>();
	FilterDataPoints(target, *path, radius, level, delete_points, results);
	delete path;
}

template<typename T,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPLeaf<T,BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
	std::vector<mvp::datapoint_t<T,PL>> results;
	for (int i=0;i < this->m_nvps;i++){
		vp_t<T> vp = this->m_vps[i];
		if (vp.active){
			results.push_back({ vp.id, vp.key });
		}
	}
	for (int i=0;i < m_npoints;i++){
		if (m_points[i].active)
			results.push_back(m_points[i]);
	}
	return results;
}



#endif /* _MVPNODE_H */
