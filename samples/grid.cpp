/*
This sample illustrates how to use the grid data structure from a
minimalistic, dependency-free C++ / STL project.
*/
#include <functional>
#include <iostream>

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
std::function<double(point p)> getx = [](point p) { return p[0]; };
std::function<double(point p)> gety = [](point p) { return p[1]; };

int main(int argc, char **argv) {
  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  trajectory t2 = {{2, 2}, {4, 3}, {0, 3}};
  trajectory q = {{1, 1}, {2, 2}, {1, 3}};
  const double distThreshold = 2.0;

  frechetrange::detail::duetschvahrenhold::Grid<
      trajectory, distance_functional_type, std::function<double(point p)>,
      std::function<double(point p)>>
      grid(distThreshold, squared_dist, getx, gety);

  grid.reserve(2); // not mandatory, but adviced in case of many inserts
  grid.insert(t1);
  grid.insert(t2);
  grid.optimize();

  auto results = grid.rangeQuery(q, distThreshold);
  cout << "Data trajectories within Frechet distance " << distThreshold << ":"
       << endl;
  for (const trajectory *t : results) {
    for (const auto &p : *t)
      cout << "( " << p[0] << ", " << p[1] << " ); ";
    cout << endl;
  }

  return 0;
}
