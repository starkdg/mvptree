#ifndef _KEY_HPP
#define _KEY_HPP

struct H64KeyObject {
	uint64_t key;
	H64KeyObject():key(0){}
	H64KeyObject(const uint64_t key):key(key){}
	H64KeyObject(const H64KeyObject &other):key(other.key){}
	const H64KeyObject& operator=(const H64KeyObject &other){
		key = other.key;
		return *this;
	}
	const double distance(const H64KeyObject &other)const{
		return __builtin_popcountll(key^other.key);
	}
};


struct VectorKeyObject {
	double key[16];
	VectorKeyObject(){};
	VectorKeyObject(const double otherkey[]){
		for (int i=0;i < 16;i++) key[i] = otherkey[i];
	}
	VectorKeyObject(const VectorKeyObject &other){
		for (int i=0;i < 16;i++) key[i] = other.key[i];
	}
	const VectorKeyObject& operator=(const VectorKeyObject &other){
		for (int i=0;i < 16;i++) key[i] = other.key[i];
		return *this;
	}
	const double distance(const VectorKeyObject &other)const {
		double sum = 0;
		for (int i=0;i < 16;i++){
			sum += pow(key[i] - other.key[i], 2.0);
		}
		return sqrt(sum/16.0);
	}
	
};


#endif
