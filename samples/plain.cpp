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
typedef std::function<double(const point &, const point &)>
    distance_functional_type;
typedef std::vector<point> trajectory;
distance_functional_type squared_dist = [](const point &p1, const point &p2) {
  // provide a squared distance functional for the point type
  return (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
};
std::function<double(const point &)> getx = [](const point &p) { return p[0]; };
std::function<double(const point &)> gety = [](const point &p) { return p[1]; };

int main(int argc, char **argv) {
  //  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  //  trajectory t2 = {{1, 1}, {2, 2}, {1, 3}};
  trajectory t1 = {{0, 0}, {10, 0}}; // straight line along X-axis
  trajectory t2 = {
      {0, 1}, {5, 10}, {10, 1}}; // triangle above. Tip is 10 higher.
  // true Frechet: 10
  frechetrange::detail::duetschvahrenhold::FrechetDistance<
      distance_functional_type, decltype(getx), decltype(gety)>
      fd(squared_dist, getx, gety);

  cout << std::fixed;
  // a range scan
  for (double d = 1; d < 15; d += 0.25)
    cout << "Reachable at " << d << ":\t"
         << (fd.isBoundedBy(t1, t2, d) ? "yes" : "no") << endl;

  // estimate the distance by interval cutting
  {
    std::pair<double, double> interval = {0, 15};
    while ((interval.second - interval.first) > 0.001) {
      double m = (interval.first + interval.second) / 2;

      if (fd.isBoundedBy(t1, t2, m)) {
        interval.second = m;
      } else {
        interval.first = m;
      }
    }
    cout << "Final Interval: [" << interval.first << "," << interval.second
         << "]" << endl;
  }

  // now bringmann baldus test code, deactivated as it is not ready.
  frechetrange::detail::bringmannbaldus::FrechetDistance<
      distance_functional_type, decltype(getx), decltype(gety)>
      fd2(squared_dist, getx, gety);

  for (double d = 1; d < 15; d += 0.25)
    cout << "Reachable2 at " << d << ":\t"
         << (fd2.is_frechet_distance_at_most(t1, t2, d) ? "yes" : "no") << endl;
  {
    std::pair<double, double> interval = {0, 15};
    while ((interval.second - interval.first) > 0.001) {
      double m = (interval.first + interval.second) / 2;

      if (fd2.is_frechet_distance_at_most(t1, t2, m)) {
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
