/*
This sample illustrates how to use the Fréchet Range Queries Library from a
minimalistic, dependency-free
C++ / STL project.

*/
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "../include/frechetrange.hpp"

using std::cout;
using std::endl;

typedef std::vector<double> point;
typedef std::function<double(point, point)> distance_functional_type;
typedef std::vector<point> trajectory;
distance_functional_type squared_dist = [](point p1, point p2) {
  // provide a squared distance functional for the point type
  return (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
};
std::function<double(point p)> getx = [](point p)->double & { return p[0]; };
std::function<double(point p)> gety = [](point p) ->double & { return p[1]; };




int main(int argc, char **argv) {
  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  trajectory t2 = {{1, 1}, {2, 2}, {1, 3}};

  frechetrange::detail::duetschvahrenhold::FrechetDistance<
      distance_functional_type, std::function<double(point p)>,
      std::function<double(point p)>>
      fd(squared_dist, getx, gety);

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


    // now bringman baldus test code, deactivated as it is not ready.
/*    frechetrange::detail::bringmanbaldus::FrechetDistance<
	trajectory, trajectory::value_type, double,
	 std::function<double(point p)>, std::function<double(point p)>,distance_functional_type> fd2(
	    getx,gety,squared_dist
	);
  for (double d = 1; d < 5; d += 0.25)
    cout << "Reachable2 at " << d << ":\t"
         << (fd2.is_frechet_distance_at_most(t1,t2,d)? "yes" : "no") << endl;
*/
    
       
  return 0;
}
