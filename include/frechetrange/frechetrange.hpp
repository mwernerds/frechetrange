#ifndef FRECHET_RANGE_HPP
#define FRECHET_RANGE_HPP

/**
*   This file collects all implementations of the three winning submissions to
* the ACM SIGSPATIAL GIS Cup 2017 in a single codebase readily usable in your
* projects.
*
*   The implementations are based on simple C++ / STL concepts such that you can
* very easily write adapters to your existing data structures. See
* example/foreign.cpp for an example with a strange trajectory class.
*
*   We encourage you to use boost::geometry for all your geometry processing, as
* it provides some nice features out of the box and is also compatible with
* custom storage for geometry.
*
*   Organization of the namespaces:
*
*     frechetrange      contains boilerplate and unifying code making all
*                       relevant aspects of the submissions accessible with a
*                       single interface.
*         class euclidean_distance_sqr
*
*         detail        contains implementation details of the three
*                       approaches well-isolated
*             bb        contains implementations by Baldus and Bringmann
*                 class frechet_distance
*                 class spatial_index
*
*             bddm      contains implementations by Buchin, Diez, Diggelen, and
*                       Meulemans
*                 class spatial_hash
*
*             dv        contains implementations by Dütsch and Vahrenhold
*                 class frechet_distance
*                 class grid
*
*/
#include "distance_sqr.hpp"

#include "detail/bb/frechet_distance.hpp"
#include "detail/bb/spatial_index.hpp"

#include "detail/bddm/spatial_hash.hpp"

#include "detail/dv/frechet_distance.hpp"
#include "detail/dv/grid.hpp"

#endif
