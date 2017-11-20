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


# Contributors
Martin Werner - boilerplate, R-package, some integration work
More to come.
