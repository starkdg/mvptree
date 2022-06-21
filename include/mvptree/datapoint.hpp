#ifndef _DATAPOINT_H
#define _DATAPOINT_H


namespace mvp {
	
	template<typename T>
	struct item_t {
		long long id;
		T key;
		item_t(){}
		item_t(const long long id, const T key):id(id),key(key){}
		item_t(const item_t<T> &other):id(other.id),key(other.key){}
		item_t<T>& operator=(const item_t<T> &other){
			id = other.id;
			key = other.key;
			return *this;
		}
		const double distance(const item_t<T> &other)const{
			return key.distance(other.key);
		}
	};


	template<typename T>
	struct vp_t {
		long long id;
		T key;
		bool active;
		vp_t(){}
		vp_t(const long long id, const T key):id(id),key(key),active(true){}
		vp_t(const vp_t<T> &other):id(other.id),key(other.key),active(other.active){}
		const vp_t<T>& operator=(const vp_t<T> &other){
			id = other.id;
			key = other.key;
			active = other.active;
			return *this;
		}
		const double distance(const T &other)const{
			return key.distance(other.key);
		}
	};


	template<typename T, int PL>
	struct datapoint_t {
		static unsigned long n_query_ops;
		static unsigned long n_build_ops;
		long long id;
		T key;
		bool active;
		double dists[PL];
		datapoint_t():id(0),active(true){
			for (int i=0;i < PL;i++) dists[i] = 0;
		}
		datapoint_t(const long long id, T key):id(id),key(key),active(true){
			for (int i=0;i < PL;i++) dists[i] = 0;
		}
		datapoint_t(const datapoint_t<T, PL> &other):id(other.id),key(other.key),active(other.active){
			for (int i=0;i < PL;i++) dists[i] = other.dists[i];
		}

		datapoint_t<T,PL>& operator=(const datapoint_t<T,PL> &other){
			id = other.id;
			key = other.key;
			active = other.active;
			for (int i=0;i < PL;i++) dists[i] = other.dists[i];
			return *this;
		}

		const double distance(const datapoint_t<T,PL> &other)const{
			return key.distance(other.key);
		}
	
		const double distance(const item_t<T> &other)const{
			return key.distance(other.key);
		}
		const double distance(const vp_t<T> &other)const{
			datapoint_t<T,PL>::n_build_ops++;
			return key.distance(other.key);
		}
		const double distance(const T &other){
			datapoint_t<T,PL>::n_query_ops++;
			return key.distance(other.key);
		}
	};
}

template<typename T, int PL>
unsigned long mvp::datapoint_t<T,PL>::n_query_ops = 0;

template<typename T, int PL>
unsigned long mvp::datapoint_t<T,PL>::n_build_ops = 0;

#endif
