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

#ifndef _KEY_HPP
#define _KEY_HPP
#include <cstdint>

struct H64KeyObject {
	std::uint64_t key;
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
