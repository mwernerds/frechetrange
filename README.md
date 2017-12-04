# frechetrange
Range Queries under True Fréchet Distance - Consolidated Outcomes of the ACM SIGSPATIAL GIS Cup 2017

# Scope
This repository aims to consolidate and bring together the implementations of the three winners of the ACM SIGSPATIAL Cup 2017 
on range queries under Fréchet distance in a directly usable form.

We present a single-file, header-only unified implementation relying on boost::geometry, some examples, as well as an R package 
that we intend to bring to CRAN, soon.

Remark: THIS IS WORK IN PROGRESS. NOTHING WORKS NOW ;-)
# Concept and Basic Usage

The implementations in this library are kept as generic as possible. We are using free functions instead of object-orientation and member
functions as much as possible and we support any type of storage as long as efficient free functions for the needed data operations are
made accessible.

While we plan a highly accessible layer simplifying syntax, we are just working on basic algorithm integration. The algorithms origin from
three different groups in completely different context. Hence, this is not trivial and will take some time. However, the current
trajectory concept is as follows:

A trajectory can be represented as an arbitrary C++ class as long as it supports the following parts of the container concept:

Given a trajectory T:
- size_t T.size(): Returns the size of the trajectory (in constant time!)
- T[i]: Returns the i-th trajectory sample (in constant time).

We plan to expect the following aspects later, but they are in fact not yet needed:
- T is copy-constructible
- T provides random access iterators (begin, end, etc.)

All other aspects of the algorithm are expected to be specified as lambda-expressions or free functions. For example,
the Fréchet decider of Dütsch and Vahrenhold uses the following methods relating points:

- A squared, fast distance functional (point, point) -> (double or similar)
- A getter for the X coordinate (point) -> (double or similar)
- A getter for the Y coordinate (point) -> (double or similar)

In addition to that, we currently expect that the geometry is Euclidean.
#Example

For now, only the details interfaces are available. You need to specify explicitly all aspects of your data type as for example in the
following line from sample/bg.cpp
```
    frechetrange::detail::duetschvahrenhold::FrechetDistance<
      std::function<double(point_type, point_type)>, // the squared distance signature
      std::function<double(point_type)>, // the X getter signature
      std::function<double(point_type)> > // the Y getter signature
      fd(
	  [](point_type p1, point_type p2) {return bg::comparable_distance(p1,p2);}, // the squared distance
	  [](point_type p){return bg::get<0>(p);},
	  [](point_type p){return bg::get<1>(p);}
      );

```
As you can see, we give lambda-expressions for the three needed functionals squared_distance (comparable_distance in boost::geometry is some
fast implementation of a distance such that compare is correct). It is the squared distance for Euclidean geometries, but might be something
different for other geometries. Additionally, we forward the boost::geometry getter for X and Y as the second and third constructor parameters.

After constructing this object, it can be used quite cleanly, for example as in
```
  // estimate the distance by interval cutting
  std::pair<double, double> interval = {0, 5};        
  while ((interval.second - interval.first) > 0.001) {
    double m = (interval.first + interval.second) / 2;

    if (fd.isBoundedBy(t1, t2, m)) {
      interval.second = m;
    } else {
      interval.first = m;
    }
  }
  cout << "Final Interval: [" << interval.first << "," << interval.second << "]"
       << endl;

```

# Contributors
Martin Werner - boilerplate, R-package, some integration work
More to come.

Fabian Dütsch - Grid, optimization of the Fréchet decider
