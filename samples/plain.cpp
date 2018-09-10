/*
This sample illustrates how to use the Fréchet Range Queries Library from a
minimalistic, dependency-free
C++ / STL project.

*/
#include <cmath>
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
    template <size_t dim>
    static double get(const point &p) {
        return std::get<dim>(p);
    }
};

int main(int argc, char **argv) {
    //  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
    //  trajectory t2 = {{1, 1}, {2, 2}, {1, 3}};
    trajectory t1 = {{0, 0}, {10, 0}};  // straight line along X-axis
    trajectory t2 = {
        {0, 1}, {5, 10}, {10, 1}};  // triangle above. Tip is 10 higher.
    // true Frechet: 10
    frechetrange::detail::dv::frechet_distance<
        2, get_coordinate,
        frechetrange::euclidean_distance_sqr<2, get_coordinate>>
        fd;

    cout << std::fixed;
    // a range scan
    for (double d = 1; d < 15; d += 0.25)
        cout << "Reachable at " << d << ":\t"
             << (fd.is_bounded_by(t1, t2, d) ? "yes" : "no") << endl;

    // estimate the distance by interval cutting
    {
        std::pair<double, double> interval = {0, 15};
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

    // now bringmann baldus test code.
    frechetrange::detail::bb::frechet_distance<
        2, get_coordinate,
        frechetrange::euclidean_distance_sqr<2, get_coordinate>>
        fd2;

    for (double d = 1; d < 15; d += 0.25)
        cout << "Reachable2 at " << d << ":\t"
             << (fd2.is_bounded_by(t1, t2, d) ? "yes" : "no") << endl;
    {
        std::pair<double, double> interval = {0, 15};
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
