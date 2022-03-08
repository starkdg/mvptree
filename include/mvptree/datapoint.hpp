#ifndef _DATAPOINT_H
#define _DATAPOINT_H

#include <cmath>

using namespace std;

class DataPoint {
protected:
	long long id;
	bool active;
public:

	static int n_ops;
	
	DataPoint(const long long id):id(id),active(true){}

	DataPoint(const DataPoint &other){
		id = other.id;
		active = other.active;
	}
	
	virtual ~DataPoint(){}
	
	DataPoint& operator=(const DataPoint &other){
		id = other.id;
		active = other.active;
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

	static void ResetCount(){
		DataPoint::n_ops = 0;
	}

	static const int GetCount(){
		return DataPoint::n_ops;
	}
	
	virtual const double distance(const DataPoint *other)const = 0;
	
};

int DataPoint::n_ops = 0;

class DblDataPoint : public DataPoint {
protected:
	double value;
public:
	DblDataPoint(const long long id, const double value):DataPoint(id),value(value){}
	DblDataPoint(const DblDataPoint &other):DataPoint(other){
		value = other.value;
	}
	~DblDataPoint(){}
	
	DblDataPoint& operator=(const DblDataPoint &other){
		DataPoint::operator=(other);
		value = other.value;
		return *this;
	}

	
	const double GetValue()const{
		return value;
	}

	const double distance(const DataPoint *other)const override{
		DataPoint::n_ops++;
		DblDataPoint *pnt = (DblDataPoint*)other;
		return abs(this->value - pnt->value);
	}
};


class H64DataPoint : public DataPoint {
protected:
	uint64_t value;
public:
	H64DataPoint(const long long id, const uint64_t value):DataPoint(id),value(value){}
	H64DataPoint(const H64DataPoint &other):DataPoint(other){
		value = other.value;
	}
	~H64DataPoint(){}
	H64DataPoint& operator=(const H64DataPoint &other){
		DataPoint::operator=(other);
		value = other.value;
		return *this;
	}

	const uint64_t GetValue(){
		return value;
	}

	const double distance(const DataPoint *other)const override{
		DataPoint::n_ops++;
		H64DataPoint *pnt = (H64DataPoint*)other;
		uint64_t tmp = pnt->GetValue();
		return __builtin_popcountll(value^tmp);
	}
};

#endif
