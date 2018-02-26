/*
This sample illustrates how to use the Fréchet Range Queries Library from a
typical Boost::Geometry project.

two versions are given:
  one is vector_of_points, which models a trajectory as a vector of boost points
  one is linestring, which models a trajectory as a linestring feature.

*/
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include "../include/frechetrange.hpp"

using std::cout;
using std::endl;
namespace bg = boost::geometry;

// Declare Geometry
typedef bg::model::point<double, 2, bg::cs::cartesian> point_type;
typedef bg::model::linestring<point_type> linestring_type;

struct get_coordinate {
  template <size_t dim> static double get(const point_type &p) {
    return p.get<dim>();
  }
};

struct squared_distance {
  double operator()(const point_type &p, const point_type &q) const {
    return bg::comparable_distance(p, q);
  }
};

int main(int argc, char **argv) {
  // Linestring version
  linestring_type t1, t2;

  // read from WKT
  bg::read_wkt("LINESTRING(0 0 , 0 1, 0 2)", t1);
  bg::read_wkt("LINESTRING(1 1, 2 2, 1 3)", t2);

  // instantiate some decider, this time fully specified.
  frechetrange::detail::duetschvahrenhold::frechet_distance<2, get_coordinate,
                                                            squared_distance>
      fd;

  cout << std::fixed;
  // a range scan
  for (double d = 1; d < 5; d += 0.25)
    cout << "Reachable at " << d << ":\t"
         << (fd.is_bounded_by(t1, t2, d) ? "yes" : "no") << endl;

  // estimate the distance by interval cutting
  {
    std::pair<double, double> interval = {0, 5};
    while ((interval.second - interval.first) > 0.001) {
      double m = (interval.first + interval.second) / 2;

      if (fd.is_bounded_by(t1, t2, m)) {
        interval.second = m;
      } else {
        interval.first = m;
      }
    }
    cout << "Final Interval: [" << interval.first << "," << interval.second
         << "]" << endl;
  }

  // Bringmann Baldus Case

  frechetrange::detail::baldusbringmann::frechet_distance<2, get_coordinate,
                                                          squared_distance>
      fd2;

  cout << std::fixed;
  // a range scan
  for (double d = 1; d < 5; d += 0.25)
    cout << "Reachable at " << d << ":\t"
         << (fd2.is_bounded_by(t1, t2, d) ? "yes" : "no") << endl;

  // estimate the distance by interval cutting
  {
    std::pair<double, double> interval = {0, 5};
    while ((interval.second - interval.first) > 0.001) {
      double m = (interval.first + interval.second) / 2;

      if (fd2.is_bounded_by(t1, t2, m)) {
        interval.second = m;
      } else {
        interval.first = m;
      }
    }
    cout << "Final Interval: [" << interval.first << "," << interval.second
         << "]" << endl;
  }

  return 0;
}
