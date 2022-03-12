#ifndef _DATAPOINT_H
#define _DATAPOINT_H

#include <cmath>
#include <cstdint>

using namespace std;

class DataPoint {
protected:
	long long id;
	uint64_t value;
	bool active;
	
public:

	static int n_ops;
	
	DataPoint(const long long id, const uint64_t value):id(id),value(value),active(true){}

	DataPoint(const DataPoint &other){
		id = other.id;
		value = other.value;
		active = other.active;
	}
	
	~DataPoint(){}
	
	DataPoint& operator=(const DataPoint &other){
		id = other.id;
		value = other.value;
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

	const uint64_t GetValue()const{
		return value;
	}
	
	const double distance(const DataPoint *other)const{
		DataPoint::n_ops++;
		return __builtin_popcountll(value^(other->value));
	}

	const double distance(const uint64_t othervalue)const{
		DataPoint::n_ops++;
		return __builtin_popcountll(value^othervalue);
	}
};

int DataPoint::n_ops = 0;


#endif
