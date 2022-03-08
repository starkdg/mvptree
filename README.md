# MVPTree (Multiple Vantage Point Tree Data structure)

A distance-based tree data structure for the indexing of point
data.  Points are indexed according to its distance from selected
vantage points.  The vantage points are arbitrarily taken from the
data itself.  The data points are simply scalar, but the code
can be easily modified to work for any multi-dimensional data in
a proper metric space.  In other words, it must have a relevant
distance function `distance(x, y) >= 0` for which the following three
axioms hold:

For any points, x, y and z,

1. distance(x, y) == 0 iff x == y  (identity)
2. distance(x, y) == distance(y, x) (symmetry)
2. distance(x, y) <= distance(x,z) + distance(z, y) (triangle inequality)

## Use 

DataPoint is an abstract class intended to hold individual data points.  Subclass
from here to hold application specific data.  Two existing DataPoint subclasses
are DblDataPoint and H64DataPoint.

DblDataPoint is a general demonstration.  It represents a single floating point
value.

H64DataPoint represents a 64-bit unsigned integer and uses the hamming distance
which exists in the discrete 64-dim hamming space.  Distance is calculated by
the number of bits that are different from each other.  


Distances are simply the absolute value of the difference between two values.
The value is just a scalar but can be easily modified to be any multi-dimensional
vector with a corresponding metric distance.  It must be a proper metric distance,
obeying the following:


The interface is largely self-explanatory.  Here are the relevant
methods for the MVPTree class:

```
void Add(DataPoint *dp);
void Add(vector<DataPoint*> &points)
```

Datapoints can be added separately, but it is better to add points in batches.  For that
reason, the single Add method adds each point to a buffer and invokes the batch method when
a threshold number is accumulated.    

```
void Sync()
```

adds the buffer to the data structure before reaching the threshold.


```
const vector<DataPoint*> Query(const DataPoint &target, const double radius)
```

queries a target DataPoint for all points that fall within a radius of the target.


```
void Delete(long long id)
```

deletes the datapoint.

```
const int Size()
```

returns the number of DataPoints currently indexed.

```
const size_t MemoryUsage()
```

returns the number of bytes the index currently occupies.

```
void Clear()
```

clears the index.

See `mvptree.hpp` header for fuller definition of the api.

## Install

```
cmake .
make
make install
```

You can also run `ctest` or `make test` to run the tests.  



