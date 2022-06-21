# MVPTree (Multiple Vantage Point Tree Data structure)

A distance-based data structure for the indexing of multi-dimensional data.
The MVPTree indexes vector objects according to its distance from select vantage points chosen directly from the indexed data.



## Performance Tests

The following tests measure the build and query times for various index sizes.
Data points are 16-dimensional real-valued vector object (uniform independantly distributed).
Build operations reflect the number of distance operations needed to build the tree.
Likewise, query operations are the average number of distance operations to query the index for a constant radius of 0.04.
Both build and query operations are normalized to the size of the tree and expressed as a percentage.
Each line in the chart is the average of 5 trial runs, where the tree structure is built up, queried for 10 clusters
of 10 points, and then taken down.  


|   N   |   MEM   |   Build opers   |   Build time   |   Query opers   |   Query time   |
|-------|---------|-----------------|----------------|-----------------|----------------|
| 100K  |  122MB  |  2554%  |  21.4$mu;s  |  0.0087%  |  42.4$mu;s  |
| 200K  |  289MB  |  2693%  |  26.3$mu;s  |  0.0048%  |  63.4&mu;s  |   
| 400K  |  548MB  |  2792%  |  36.8$mu;s  |  0.0022%  |  70.2&mu;s  |
| 800K  |  963MB  |  2872%  |  59.4&mu;s  |  0.0012%  |  86.0&mu;s  |
|  1M   |  1.1GB  |  2899%  |  62.9&mu;s  |  0.0010%  |  107&mu;s |
|  2M   |  2.7GB  |  3021%  |  79.6&mu;s  |  0.0005%  |  163&mu;s |
|  4M   |  5.5GB  |  3139%  |  101&mu;s   |  0.0002%  |  242&mu;s |



Here are the query stats for various radius sizes and a constant index size of N = 1,000,000:

| radius | Query Opers. | Query Time |
|--------|--------------|------------|
| 0.02 | 0.00094% | 30.9&mu;s |
| 0.04 | 0.00092% | 121&mu;s |
| 0.06 | 0.0012% | 359&mu;s |
| 0.08 | 0.00041% | 971&mu;s |
| 0.10 | 0.0241% | 2.67&mu;s |


To try it for yourself, run the test program, `runmvptree2`.  You can fiddle with the test parameters
in the file, tests/run_mvptree2.cpp.  

## Programming

Use the template header files.
First, create a custom data type. Examples are provided under key.hpp:

```
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
```

The custom data type must provide the copy contructor, assignment operator and distance function implementation.  


Next, use the templated class:

```
#include "mvptree/mvptree.hpp"


MVPTree<VectorKeyObject> mvptree;

vector<VectorKeyObject> points;

mvptree.Add(points);
mvptree.Sync();

int sz = mvptree.Size();
assert(sz == points.size());


VectorKeyObject target;
double radius = ...  //  value for data set 
vector<item_t<VectorKeyObject>> results = mvptree.Query(target, radius);


mvptree.Clear();
sz = mvptree.Size();
assert(sz == 0);

```

For more details, see the programs in the tests directory. 

## Install

```
cmake .
make
make test
make install
```





