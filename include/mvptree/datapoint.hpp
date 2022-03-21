#ifndef _DATAPOINT_H
#define _DATAPOINT_H

#include <cstring>
#include <stdexcept>
#include <cmath>

using namespace std;

template<int PL>
class DataPoint {
protected:
	double dists[PL];
	long long id;
	bool active;
public:

	static int n_ops;
	
	DataPoint(const long long id):id(id),active(true){}

	DataPoint(const DataPoint<PL> &other){
		id = other.id;
		active = other.active;
		memcpy(dists, other.dists, PL*sizeof(double));
	}
	
	virtual ~DataPoint(){}
	
	DataPoint<PL>& operator=(const DataPoint<PL> &other){
		id = other.id;
		active = other.active;
		memcpy(dists, other.dists, PL*sizeof(double));
		return *this;
	}

	const bool IsActive() const{
		return active;
	}
	
	void Deactivate(){
		active = false;
	}

	const long long GetId()const{
		return id;
	}

	const double GetPath(const int n)const{
		if (n < 0 || n >= PL) throw out_of_range("index out of range");
		return dists[n];
	}

	void SetPath(const double val, const int n){
		if (n < 0 || n >= PL) throw out_of_range("index out of range");
		dists[n] = val;
	}
	
	virtual const double distance(const DataPoint<PL> *other)const = 0;
	
};

template<int PL>
class DblDataPoint : public DataPoint<PL> {
protected:
	double value;
public:
	DblDataPoint(const long long id, const double value):DataPoint<PL>(id),value(value){
		this->value = value;
	}
	DblDataPoint(const DblDataPoint &other):DataPoint<PL>(other){
		value = other.value;
	}
	~DblDataPoint(){}
	
	DblDataPoint& operator=(const DblDataPoint &other){
		DataPoint<PL>::operator=(other);
		value = other.value;
		return *this;
	}

	
	const double GetValue()const{
		return value;
	}

	const double distance(const DataPoint<PL> *other)const override{
		DataPoint<PL>::n_ops++;
		DblDataPoint *pnt = (DblDataPoint*)other;
		return abs(this->value - pnt->value);
	}
};

template<int PL>
class H64DataPoint : public DataPoint<PL> {
protected:
	uint64_t value;
public:
	H64DataPoint(const long long id, const uint64_t value):DataPoint<PL>(id),value(value){}
	H64DataPoint(const H64DataPoint &other):DataPoint<PL>(other){
		value = other.value;
	}
	~H64DataPoint(){}
	H64DataPoint& operator=(const H64DataPoint &other){
		DataPoint<PL>::operator=(other);
		value = other.value;
		return *this;
	}

	const uint64_t GetValue()const{
		return value;
	}

	const double distance(const DataPoint<PL> *other)const override{
		DataPoint<PL>::n_ops++;
		H64DataPoint *pnt = (H64DataPoint*)other;
		uint64_t tmp = pnt->GetValue();
		return __builtin_popcountll(value^tmp);
	}
};


template<int PL>
int DataPoint<PL>::n_ops = 0;

#endif
