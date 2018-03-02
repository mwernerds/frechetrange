/*
This sample illustrates how to use the grid data structure from a
minimalistic, dependency-free C++ / STL project.
*/
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "../include/frechetrange/frechetrange.hpp"

using std::cout;
using std::endl;

typedef std::pair<double, double> point;
typedef std::vector<point> trajectory;

struct get_coordinate {
  template <size_t dim> static double get(const point &p) {
    return std::get<dim>(p);
  }
};

int main(int argc, char **argv) {
  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  trajectory t2 = {{2, 2}, {4, 3}, {0, 3}};
  trajectory q1 = {{1.5, 1.5}, {3, 2.5}, {3.5, 3.5}, {0, 2}};
  const double distThreshold1 = 1.5;
  trajectory q2 = {{1, 1}, {2, 2}, {1, 3}};
  const double distThreshold2 = 2.0;

  frechetrange::detail::bddm::spatial_hash<2, trajectory, get_coordinate>
      spatialHash;

  spatialHash.insert(t1);
  spatialHash.insert(t2);
  spatialHash.build_index();

  // first version of rangeQuery: returning the result set
  auto results = spatialHash.range_query(q1, distThreshold1);

  cout << "Data trajectories within Frechet distance " << distThreshold1
       << " of q1:" << endl;
  for (const trajectory *t : results) {
    for (const point &p : *t)
      cout << "( " << p.first << ", " << p.second << " ); ";
    cout << endl;
  }

  // second version of rangeQuery: using an output functional
  cout << endl;
  cout << "Data trajectories within Frechet distance " << distThreshold2
       << " of q2:" << endl;
  auto output = [](const trajectory &t) {
    for (const point &p : t)
      cout << "( " << p.first << ", " << p.second << " ); ";
    cout << endl;
  };
  spatialHash.range_query(q2, distThreshold2, output);

  return 0;
}
