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
typedef std::function<double(const point &, const point &)>
    distance_functional_type;
typedef std::vector<point> trajectory;

distance_functional_type squared_dist = [](point p1, point p2) {
  // provide a squared distance functional for the point type
  return (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
};
std::function<double(const point &)> getx = [](const point &p) { return p[0]; };
std::function<double(const point &)> gety = [](const point &p) { return p[1]; };

int main(int argc, char **argv) {
  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  trajectory t2 = {{2, 2}, {4, 3}, {0, 3}};
  trajectory q1 = {{1.5, 1.5}, {3, 2.5}, {3.5, 3.5}, {0, 2}};
  const double distThreshold1 = 1.5;
  trajectory q2 = {{1, 1}, {2, 2}, {1, 3}};
  const double distThreshold2 = 2.0;

  double meshSize = std::max(distThreshold1, distThreshold2);
  frechetrange::detail::duetschvahrenhold::Grid<
      trajectory, distance_functional_type, decltype(getx), decltype(gety)>
      grid(meshSize, squared_dist, getx, gety);

  grid.reserve(2); // not mandatory, but advised in case of many inserts
  grid.insert(t1);
  grid.insert(t2);
  grid.optimize(); // not mandatory, but advised after completing inserts

  // first version of rangeQuery: returning the result set
  auto results = grid.rangeQuery(q1, distThreshold1);

  cout << "Data trajectories within Frechet distance " << distThreshold1
       << " of q1:" << endl;
  for (const trajectory *t : results) {
    for (const point &p : *t)
      cout << "( " << p[0] << ", " << p[1] << " ); ";
    cout << endl;
  }

  // second version of rangeQuery: using an output functional
  cout << endl;
  cout << "Data trajectories within Frechet distance " << distThreshold2
       << " of q2:" << endl;
  auto output = [](const trajectory &t) {
    for (const point &p : t)
      cout << "( " << p[0] << ", " << p[1] << " ); ";
    cout << endl;
  };
  grid.rangeQuery(q2, distThreshold2, output);

  return 0;
}
