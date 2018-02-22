/*
This sample illustrates how to use the grid data structure from a
minimalistic, dependency-free C++ / STL project.
*/
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "../include/tue.hpp"

using std::cout;
using std::endl;

typedef std::pair<double, double> point;
typedef std::vector<point> trajectory;

std::function<double(const point &)> getx = [](const point &p) {
  return std::get<0>(p);
};
std::function<double(const point &)> gety = [](const point &p) {
  return std::get<1>(p);
};

int main(int argc, char **argv) {
  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
  trajectory t2 = {{2, 2}, {4, 3}, {0, 3}};
  trajectory q1 = {{1.5, 1.5}, {3, 2.5}, {3.5, 3.5}, {0, 2}};
  const double distThreshold1 = 1.5;
  trajectory q2 = {{1, 1}, {2, 2}, {1, 3}};
  const double distThreshold2 = 2.0;

  tue_details::SpatialHash<trajectory, decltype(getx), decltype(gety)>
      spatialHash(getx, gety);

  spatialHash.insert(t1);
  spatialHash.insert(t2);
  spatialHash.buildIndex();

  // first version of rangeQuery: returning the result set
  auto results = spatialHash.rangeQuery(q1, distThreshold1);

  cout << "Data trajectories within Frechet distance " << distThreshold1
       << " of q1:" << endl;
  for (const trajectory *t : results) {
    for (const point &p : *t)
      cout << "( " << getx(p) << ", " << gety(p) << " ); ";
    cout << endl;
  }

  // second version of rangeQuery: using an output functional
  cout << endl;
  cout << "Data trajectories within Frechet distance " << distThreshold2
       << " of q2:" << endl;
  auto output = [](const trajectory &t) {
    for (const point &p : t)
      cout << "( " << getx(p) << ", " << gety(p) << " ); ";
    cout << endl;
  };
  spatialHash.rangeQuery(q2, distThreshold2, output);

  return 0;
}
