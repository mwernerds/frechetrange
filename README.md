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
- T is copy-constructible
- size_t T.size(): Returns the size of the trajectory (in constant time!)
- T[i]: Returns the i-th trajectory sample (in constant time).

All other aspects of the algorithm are expected to be specified as template parameters and functionals. For example,
the Fréchet deciders require the following parameters:

- The number of dimensions of the points
- A getter for each dimension's coordinate get<dim>(const point&) -> (double or similar)
- A squared, fast distance functional (const point&, const point&) -> (double or similar)

In addition to that, we currently expect that the geometry is Euclidean.
#Example

For now, only the details interfaces are available. You need to at least specify the number of dimensions and the getter type of your
point type explicitly as in the following example (compare to sample/bg.cpp):
```
    frechetrange::detail::dv::frechet_distance<
      2,											         	   // number of dimensions
      get_coordinate,                                        	   // such that get_coordinate::get<dim>(p) returns p's dim-th coordinate
      std::function<double(const point_type&, const point_type&)>> // the squared distance signature
      fd(
	  [](const point_type& p, const point_type& q) {
	    return bg::comparable_distance(p1,p2);                     // the squared distance functional
      });

```
As you can see, we can give a lambda-expression for the squared distance functional (comparable_distance in boost::geometry is some
fast implementation of a distance such that compare is correct). It is the squared distance for Euclidean geometries, but might be something
different for other geometries, in later versions.

After constructing this object, it can be used quite cleanly, for example as in
```
  // estimate the distance by interval cutting
  std::pair<double, double> interval = {0, 5};        
  while ((interval.second - interval.first) > 0.001) {
    double m = (interval.first + interval.second) / 2;

    if (fd.is_bounded_by(t1, t2, m)) {
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

Fabian Dütsch - adapting the implementations, optimizations
