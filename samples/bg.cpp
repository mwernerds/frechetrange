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




int main(int argc, char **argv) {

    // Linestring version
    
    linestring_type t1,t2;

    // read from WKT
    bg::read_wkt("LINESTRING(0 0 , 0 1, 0 2)",t1);     
    bg::read_wkt("LINESTRING(1 1, 2 2, 1 3)",t2);     

    // instantiate some decider, this time fully specified.
    frechetrange::detail::duetschvahrenhold::FrechetDistance<
      std::function<double(point_type, point_type)>, // the squared distance signature
      std::function<double(point_type)>, // the X getter signature
      std::function<double(point_type)> > // the Y getter signature
      fd(
	  [](point_type p1, point_type p2) {return bg::comparable_distance(p1,p2);}, // the squared distance
	  [](point_type p){return bg::get<0>(p);},
	  [](point_type p){return bg::get<1>(p);}
      );

  cout << std::fixed;
  // a range scan
  for (double d = 1; d < 5; d += 0.25)
    cout << "Reachable at " << d << ":\t"
         << (fd.isBoundedBy(t1, t2, d) ? "yes" : "no") << endl;

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

  return 0;
}
