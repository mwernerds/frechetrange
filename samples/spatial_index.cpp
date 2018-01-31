/*
This sample illustrates how to use the spatial_index data structure from a
minimalistic, dependency-free C++ / STL project.
*/
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

  frechetrange::detail::bringmannbaldus::spatial_index<
      trajectory, distance_functional_type, decltype(getx), decltype(gety)>
      spatial_index(squared_dist, getx, gety);

  spatial_index.add_curve(1, t1);
  spatial_index.add_curve(2, t2);

  std::vector<int> results1 =
      spatial_index.get_close_curves(q1, distThreshold1);
  cout << "Data trajectories within Frechet distance " << distThreshold1
       << " of q1:" << endl;
  for (int i : results1) {
    cout << "t" << i << endl;
  }

  std::vector<int> results2 =
      spatial_index.get_close_curves(q2, distThreshold2);
  cout << "Data trajectories within Frechet distance " << distThreshold2
       << " of q2:" << endl;
  for (int i : results2) {
    cout << "t" << i << endl;
  }

  return 0;
}
