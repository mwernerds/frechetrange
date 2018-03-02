#ifndef FRECHET_RANGE_HPP
#define FRECHET_RANGE_HPP

/**
*   This file collects all implementations of the three winning submissions to
* the ACM SIGSPATIAL GIS Cup 2017 in a single codebase readily usable in your
* projects.
*
*   The implementations are based on simple C++ / STL concepts such that you can
* very easily write adapters to your existing data structures. See
* example/foreign.cpp for an example with a strange trajectory class.
*
*   We encourage you to use boost::geometry for all your geometry processing, as
* it provides some nice features out of the box and is also compatible with
* custom storage for geometry.
*
*   Organization of this file:
*
*   Namespaces:
*
*     frechetrange      contains boilerplate and unifying code making all
*                       relevant aspects of the submissions accessible with a
*                       single interface.
*         struct euclidean_distance_sqr
*         detail            contains implementation details of the three
*                           approaches well-isolated
*             duetschvarhenhold
*                 class frechet_distance
*                 class grid
*             baldusbringmann
*                 class frechet_distance
*                 class spatial_index
*
*
*/

#include <algorithm> // for std::sort, std::lower_bound, std::upper_bound, and std::max
#include <array>
#include <cmath>      // for std::floor
#include <cstddef>    // for std::size_t
#include <functional> // for std::hash
#include <memory>     // for std::unique_ptr
#include <stdexcept>  // for std::invalid_argument
#include <unordered_map>
#include <utility> // for std::move
#include <vector>

#ifdef ENABLE_MULTITHREADING
#include <deque>
#include <future> // for std::async
#endif

namespace frechetrange {
namespace detail {
// template meta programming
template <typename point, typename get_coordinate, size_t dim>
struct sum_of_sqr_diffs {
  static double dist(const point &p, const point &q) {
    auto delta_dim = get_coordinate::template get<dim - 1>(p) -
                     get_coordinate::template get<dim - 1>(q);
    return static_cast<double>(delta_dim * delta_dim) +
           sum_of_sqr_diffs<point, get_coordinate, dim - 1>::dist(p, q);
  }
};

template <typename point, typename get_coordinate>
struct sum_of_sqr_diffs<point, get_coordinate, 0> {
  static double dist(const point &, const point &) { return 0.0; }
};

template <typename point, typename get_coordinate, size_t dim>
struct add_sum_of_pairwise_diff_prods {
  static void calc(const point &p, const point &q, const point &m, double &a,
                   double &b, double &c) {
    // u := q - p
    auto u_dim = get_coordinate::template get<dim - 1>(q) -
                 get_coordinate::template get<dim - 1>(p);
    // v := p - m
    auto v_dim = get_coordinate::template get<dim - 1>(p) -
                 get_coordinate::template get<dim - 1>(m);

    a += u_dim * u_dim;
    b += u_dim * v_dim;
    c += v_dim * v_dim;

    add_sum_of_pairwise_diff_prods<point, get_coordinate, dim - 1>::calc(
        p, q, m, a, b, c);
  }
};

template <typename point, typename get_coordinate>
struct add_sum_of_pairwise_diff_prods<point, get_coordinate, 0> {
  static void calc(const point &, const point &, const point &, double &a,
                   double &b, double &c) {}
};
}

template <size_t dimensions, typename get_coordinate>
struct euclidean_distance_sqr {
  template <typename point>
  double operator()(const point &p, const point &q) const {
    return detail::sum_of_sqr_diffs<point, get_coordinate, dimensions>::dist(p,
                                                                             q);
  }
};

namespace detail {
namespace duetschvahrenhold {
// source code from Fabian Dütsch, templatized and generalized by Martin
/**
* Arbitrary constant > 1.0 to signal that a segment is not reachable.
*/
static constexpr double BEGIN_NOT_REACHABLE = 2.0;

template <size_t dimensions, typename get_coordinate,
          typename squared_distance =
              euclidean_distance_sqr<dimensions, get_coordinate>>
class frechet_distance {
public:
  frechet_distance(const squared_distance &dist2 = squared_distance())
      : _dist2(dist2), _leftSegmentBegins() {}

  frechet_distance(const frechet_distance &) = default;
  frechet_distance(frechet_distance &&) = default;
  frechet_distance &operator=(const frechet_distance &) = default;
  frechet_distance &operator=(frechet_distance &&) = default;

  /**
  * Returns whether the Fréchet distance of the passed trajectories is bounded
  * by the passed upper bound.
  * Not thread-safe.
  * @pre Neither of the trajectories is empty.
  */
  template <typename Trajectory>
  bool is_bounded_by(const Trajectory &traj1, const Trajectory &traj2,
                     double distance_bound) const {
    const double boundSquared = distance_bound * distance_bound;

    // ensure that the corners of the free space diagram are reachable
    if (_dist2(traj1[0], traj2[0]) > boundSquared ||
        _dist2(traj1[traj1.size() - 1], traj2[traj2.size() - 1]) > boundSquared)
      return false;

    bool firstIsSmaller = traj1.size() <= traj2.size();
    const Trajectory &smallerTraj = firstIsSmaller ? traj1 : traj2;
    const Trajectory &biggerTraj = firstIsSmaller ? traj2 : traj1;

    if (smallerTraj.size() == 1) {
      return comparePointToTrajectory(smallerTraj[0], biggerTraj, boundSquared);
    } else if (!matchInnerPointsMonotonously(traj1, traj2, boundSquared) ||
               !matchInnerPointsMonotonously(traj2, traj1, boundSquared)) {
      // There exists no monotone matchting from one trajectory to the other,
      // such that each point is matchted to a point on a segment within
      // distance.
      return false;
    }

    // use the smaller trajectory as first parameter, to consume less memory
    return traverseFreeSpaceDiagram(smallerTraj, biggerTraj, boundSquared);
  }

private:
  squared_distance _dist2;

  /**
  * Beginnings of reachable parts of the free space segments on
  * the "frontline", i. e., the right segments of the free space cells lastly
  * processed (and accordingly the left segments of the current column).
  * A segment is reachable, if its respective beginning is <= 1.0.
  */
  mutable std::vector<double> _leftSegmentBegins;

  /**
  * Decides the Fréchet distance problem for a trajectory consisting of only
  * one point.
  */
  template <typename Trajectory, typename point_type>
  bool comparePointToTrajectory(const point_type &p,
                                const Trajectory &trajectory,
                                const double boundSquared) const {
    // the first point has already been tested
    for (size_t i = 1; i < trajectory.size(); ++i) {
      if (_dist2(p, trajectory[i]) > boundSquared) {
        return false;
      }
    }
    return true;
  }

  /**
  * Ensures that the sequence of passed points can be matched to a sequence of
  * points on the passed segments,
  * such that each matching is within the passed distance bound and the
  * sequence of segment points is monotone.
  */
  template <typename Trajectory>
  bool matchInnerPointsMonotonously(const Trajectory &points,
                                    const Trajectory &segments,
                                    double boundSquared) const {
    // the last point has already been tested
    const size_t pointsToMatchEnd = points.size() - 1;
    const size_t numSegments = segments.size() - 1;

    if (pointsToMatchEnd <= 1 || numSegments == 0) {
      // nothing to test
      return true;
    }

    // the first point has already been tested
    size_t pointIdx = 1;
    size_t segIdx = 0;
    double segmentPart = 0.0;
    double begin = 0.0, end = 1.0;
    while (true) {
      // search a point on a segment within distance of the current point
      getLineCircleIntersections(segments[segIdx], segments[segIdx + 1],
                                 points[pointIdx], boundSquared, begin, end);

      if (begin <= 1.0 && end >= segmentPart) {
        // found a matching point on the current segment
        if (segmentPart < begin) {
          segmentPart = begin;
        }

        // advance to the next point
        ++pointIdx;
        if (pointIdx == pointsToMatchEnd) {
          // the last point has already been tested
          return true;
        }
      } else {
        // advance to the next segment and never go back,
        // so that the matching is monotone
        ++segIdx;
        segmentPart = 0.0;

        if (segIdx == numSegments) {
          // The end of the last segment has been reached.
          // => The remaining points cannot be matched monotonously.
          return false;
        }
      }
    }
  }

  /**
  * Returns whether the Fréchet distance of the passed trajectories is bounded
  * by the upper bound.
  * @pre Both trajectories consist of at least two points
  *      and the starting and ending points are within distance.
  */
  template <typename Trajectory>
  bool traverseFreeSpaceDiagram(const Trajectory &p1, const Trajectory &p2,
                                const double boundSquared) const {
    // init the beginnings of reachable parts of the "frontline",
    // i.e., the right segments of the column lastly processed
    // (and accordingly the left segments of the current column)
    clearFrontline(p1.size() - 1);

    // the bottom-left free space corner has been ensured to be reachable
    _leftSegmentBegins[0] = 0.0;
    // beginning of the reachable part of the bottommost segment of the current
    // column described as a scalar value
    double bottommostSegmentBegin = 0.0;
    // bottommost reachable segment of the current column
    size_t bottommostReachableRow = 0;

    const size_t lastRow = p1.size() - 2;
    const size_t lastColumn = p2.size() - 2;
    // compute the reachable parts of the free space diagram's cells columnwise
    for (size_t colIdx = 0; colIdx < lastColumn; ++colIdx) {
      // compute whether the current column's bottommost segment is reachable
      if (bottommostSegmentBegin == 0.0 &&
          _dist2(p1[0], p2[colIdx]) > boundSquared) {
        bottommostSegmentBegin = BEGIN_NOT_REACHABLE;
      }

      // beginning of the reachable part of the current bottom segment
      double currBottomSegBegin = bottommostSegmentBegin;
      double currBottomSegEnd = currBottomSegBegin;

      double currRightSegBegin, currRightSegEnd;
      // traverse this columns' reachable cells
      for (size_t rowIdx = bottommostReachableRow; rowIdx < lastRow; ++rowIdx) {
        // ensure that this cell is reachable
        if (isSegmentReachable(_leftSegmentBegins[rowIdx]) ||
            isSegmentReachable(currBottomSegBegin)) {
          // whether the cell's top right corner lies in free space
          bool isTopRightFree =
              _dist2(p1[rowIdx + 1], p2[colIdx + 1]) <= boundSquared;

          // compute the reachable part of the current cell's right
          // segment
          if (isTopRightFree && currBottomSegEnd >= 1.0 &&
              currBottomSegBegin <= 1.0) {
            // the entire right segment is reachable
            currRightSegBegin = 0.0;
            currRightSegEnd = 1.0;
          } else {
            getLineCircleIntersections(p1[rowIdx], p1[rowIdx + 1],
                                       p2[colIdx + 1], boundSquared,
                                       currRightSegBegin, currRightSegEnd);
            currRightSegBegin = getReachableBegin(
                currRightSegBegin, currRightSegEnd, _leftSegmentBegins[rowIdx],
                currBottomSegBegin);
          }

          // compute the reachable part of the current cell's top segment
          double currTopSegBegin, currTopSegEnd;
          if (isTopRightFree && _leftSegmentBegins[rowIdx + 1] <= 0.0) {
            // the entire top segment is reachable
            currTopSegBegin = 0.0;
            currTopSegEnd = 1.0;
          } else {
            getLineCircleIntersections(p2[colIdx], p2[colIdx + 1],
                                       p1[rowIdx + 1], boundSquared,
                                       currTopSegBegin, currTopSegEnd);
            currTopSegBegin = getReachableBegin(currTopSegBegin, currTopSegEnd,
                                                currBottomSegBegin,
                                                _leftSegmentBegins[rowIdx]);
          }

          // add the current cell to the "frontline"
          _leftSegmentBegins[rowIdx] = currRightSegBegin;
          currBottomSegBegin = currTopSegBegin;
          currBottomSegEnd = currTopSegEnd;
        }

        // update the bottommost reachable segment, if necessary
        if (bottommostReachableRow == rowIdx &&
            !isSegmentReachable(_leftSegmentBegins[rowIdx])) {
          ++bottommostReachableRow;
        }
      }

      // compute the reachable part of the topmost cell's right segment
      if (isSegmentReachable(_leftSegmentBegins[lastRow]) ||
          isSegmentReachable(currBottomSegBegin)) {
        if (_dist2(p1[lastRow + 1], p2[colIdx + 1]) <= boundSquared &&
            currBottomSegEnd >= 1.0 && currBottomSegBegin <= 1.0) {
          // the entire segment is reachable
          currRightSegBegin = 0.0;
        } else {
          getLineCircleIntersections(p1[lastRow], p1[lastRow + 1],
                                     p2[colIdx + 1], boundSquared,
                                     currRightSegBegin, currRightSegEnd);
          currRightSegBegin = getReachableBegin(
              currRightSegBegin, currRightSegEnd, _leftSegmentBegins[lastRow],
              currBottomSegBegin);
        }
        _leftSegmentBegins[lastRow] = currRightSegBegin;
      }

      // ensure that a segment of the frontline is reachable
      if (bottommostReachableRow == lastRow &&
          !isSegmentReachable(_leftSegmentBegins[lastRow])) {
        return false;
      }
    }

    if (isSegmentReachable(_leftSegmentBegins[lastRow])) {
      // the top-right corner is reachable via its left segment
      return true;
    }

    // traverse the rightmost column
    double currBottomSegBegin = BEGIN_NOT_REACHABLE;
    for (size_t rowIdx = bottommostReachableRow; rowIdx < lastRow; ++rowIdx) {
      // ensure that this cell is reachable
      if (isSegmentReachable(_leftSegmentBegins[rowIdx]) ||
          isSegmentReachable(currBottomSegBegin)) {
        // compute the reachable part of the current cell's top segment
        double currTopSegBegin, currTopSegEnd;
        if (_dist2(p1[rowIdx + 1], p2[lastColumn + 1]) <= boundSquared &&
            _leftSegmentBegins[rowIdx + 1] <= 0.0) {
          // the entire top segment is reachable
          currBottomSegBegin = 0.0;
        } else {
          getLineCircleIntersections(p2[lastColumn], p2[lastColumn + 1],
                                     p1[rowIdx + 1], boundSquared,
                                     currTopSegBegin, currTopSegEnd);
          currBottomSegBegin =
              getReachableBegin(currTopSegBegin, currTopSegEnd,
                                currBottomSegBegin, _leftSegmentBegins[rowIdx]);
        }
      }
    }

    // return whether the top-right corner is reachable via its bottom segment
    return isSegmentReachable(currBottomSegBegin);
  }

  /**
  * Increases the allocated memory of the frontline, if necessary,
  * and marks its segments as not reachable.
  */
  void clearFrontline(size_t rows) const {
    // reserve memory, if necessary
    if (_leftSegmentBegins.size() < rows) {
      _leftSegmentBegins.resize(rows, BEGIN_NOT_REACHABLE);
    }

    // mark the frontline as not reachable
    std::fill_n(_leftSegmentBegins.begin(), rows, BEGIN_NOT_REACHABLE);
  }

  /**
  * Computes the two intersections of the line through p1 and p2
  * with the circle around the center c with the passed radius.
  * @param[out] begin Scalar value, such that (p1 + begin * (p2-p1)) is the
  *                   first intersection, if it exists,
  *                   or BEGIN_NOT_REACHABLE, otherwise.
  * @param[out] end Scalar value, such that begin <= end and
  *                 (p1 + end * (p2-p1)) is the second intersection,
  *                 if it exists, or unchanged, otherwise.
  */
  template <typename point>
  void getLineCircleIntersections(const point &p1, const point &p2,
                                  const point &cp, const double radiusSquared,
                                  double &begin, double &end) const {
    // Compute the points p1+x*(p2-p1) on the segment from p1 to p2
    // that intersect the circle around cp with radius r:
    // d(cp, p1+x*(p2-p1))^2 = r^2  <=>  a*x^2 + bx + c = 0,
    // where a, b, c, and the auxiliary vectors u and v are defined as follows:
    // u := p2 - p1, and v := p1 - cp
    // a := u^2
    double a = 0.0;
    // b := 2 * dotProduct(u, v),
    double b = 0.0;
    // c := v^2 - r^2,
    double c = -radiusSquared;
    add_sum_of_pairwise_diff_prods<point, get_coordinate, dimensions>::calc(
        p1, p2, cp, a, b, c);
    b *= 2.0;

    // Solve by quadratic formula:
    // x_1,2 = (-b +- sqrt(b^2-4ac)) / 2a
    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
      // The circle does not intersect the line through p1 and p2.
      begin = BEGIN_NOT_REACHABLE;
      return;
    } else if (a == 0.0) {
      // => p1 = p2
      if (_dist2(p1, cp) <= radiusSquared) {
        begin = 0.0;
        end = 1.0;
      } else {
        begin = end = BEGIN_NOT_REACHABLE;
      }
      return;
    }

    double sqrtD = std::sqrt(discriminant);
    begin = (-b - sqrtD) / (2.0 * a);
    end = (-b + sqrtD) / (2.0 * a);
  }

  /**
  * Computes beginning of the reachable part of a free space segment.
  * @param currSegmentBegin The beginning of the intersection of free space
  *                         with the segment described as a scalar value.
  * @param currSegmentEnd The ending of the intersection of free space
  *                       with the segment described as a scalar value.
  * @param prevParallelSegBegin The beginning of the reachable part of the
  *                             parallel segment of the same free space cell
  *                             described as a scalar value.
  * @param prevOrthogonalSegBegin The beginning of the reachable part of the
  *                               orthogonal segment of the same free space
  *                               cell described as a scalar value.
  * @return A scalar value describing the beginning of the reachable part of
  *         the segment.
  */
  double getReachableBegin(double currSegmentBegin, double currSegmentEnd,
                           double prevParallelSegBegin,
                           double prevOrthogonalSegBegin) const {
    if (currSegmentEnd < 0.0) {
      // The segment does not intersect the free space.
      return BEGIN_NOT_REACHABLE;
    } else if (currSegmentBegin > 1.0) {
      // The segment does not intersect free space,
      // and there is no need to change currSegmentBegin.
    } else if (isSegmentReachable(prevOrthogonalSegBegin)) {
      // The reachable part equals the intersection with the free space.
    } else if (currSegmentBegin < prevParallelSegBegin) {
      // The beginning of the reachable part is greater
      // than the beginning of the intersection.
      return (currSegmentEnd >= prevParallelSegBegin) ? prevParallelSegBegin
                                                      : BEGIN_NOT_REACHABLE;
    }

    return currSegmentBegin;
  }

  /**
  * Returns whether a free space segment with the passed
  * beginning of the reachable part (as returned by getReachableBegin)
  * is reachable.
  * @param begin The beginning of the reachable part of the segment
  *              described as a scalar value.
  */
  bool isSegmentReachable(double begin) const { return begin <= 1.0; }
};

/*
* class Grid
*/
#include "duetsch_vahrenhold_grid.hpp"

} // duetschvahrenhold

#include <cassert>
#include <limits>

namespace baldusbringmann {

typedef double distance_t;

template <typename valtype> valtype sqr(valtype a) { return a * a; }

/*
 * Represents a trajectory. Additionally to the points given in the input file,
 * we also store the length of any prefix of the trajectory.
 */
template <typename Trajectory> class curve {
  const Trajectory _trajectory;
  std::vector<distance_t> _prefix_length;

public:
  template <typename squared_distance>
  curve(const Trajectory &t, const squared_distance &dist2)
      : _trajectory(t), _prefix_length(t.size()) {
    _prefix_length[0] = 0;
    for (size_t i = 1; i < t.size(); ++i)
      _prefix_length[i] =
          _prefix_length[i - 1] + std::sqrt(dist2(t[i - 1], t[i]));
  }

  curve() = default;

  size_t size() const { return _trajectory.size(); }

  const Trajectory &trajectory() const { return _trajectory; }

  distance_t curve_length(size_t i, size_t j) const {
    return _prefix_length[j] - _prefix_length[i];
  }
};

template <size_t dimensions, typename get_coordinate,
          typename squared_distance =
              euclidean_distance_sqr<dimensions, get_coordinate>>
class frechet_distance {
  static constexpr distance_t eps = 10e-10;

  squared_distance _dist2;

  /*Section 1: special types and their elementary operations*/
  typedef std::pair<distance_t, distance_t> interval; // .first is the
                                                      // startpoint, .second the
                                                      // endpoint (start
                                                      // inclusive, end
                                                      // exclusive)

  const interval empty_interval{std::numeric_limits<distance_t>::max(),
                                std::numeric_limits<distance_t>::lowest()};

  bool is_empty_interval(const interval &i) const {
    return i.first >= i.second;
  }

  template <typename point>
  interval intersection_interval(const point &circle_center, distance_t radius,
                                 const point &line_start,
                                 const point &line_end) const {
    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 =
    // radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2)
    // + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) =
    // radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 *
    // line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, b := 2 * line_start.x * v.x + 2 * line_start.y *
    // v.y, c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * b + c) = 0
    // <=> lambda^2 + (b / a) * lambda + c / a) = 0
    // <=> lambda1/2 = - (b / 2a) +/- sqrt((b / 2a)^2 - c / a)

    distance_t a = 0.0;
    distance_t b = 0.0;
    distance_t c = -sqr(radius);
    add_sum_of_pairwise_diff_prods<point, get_coordinate, dimensions>::calc(
        line_start, line_end, circle_center, a, b, c);

    if (a == 0.0) {
      // => line_start = line_end
      if (_dist2(line_start, circle_center) <= radius * radius) {
        return {0.0, 1.0};
      } else {
        return empty_interval;
      }
    }

    distance_t discriminant = sqr(b / a) - c / a;

    if (discriminant < 0) {
      return empty_interval; // no intersection;
    }

    distance_t lambda1 = -b / a - std::sqrt(discriminant);
    distance_t lambda2 = -b / a + std::sqrt(discriminant);

    if (lambda2 < 0 || lambda1 > 1)
      return empty_interval;
    else
      return {std::max<distance_t>(lambda1, 0),
              std::min<distance_t>(lambda2, 1)};
  }

  /*Section 2: Elementary Frechet operations*/
  template <typename Trajectory, typename point>
  distance_t get_dist_to_point_sqr(const Trajectory &t, const point &p) const {
    distance_t result = 0;
    for (size_t i = 0; i < t.size(); ++i)
      result = std::max(result, _dist2(t[i], p));
    return result;
  }

  template <typename Trajectory>
  interval get_reachable_a(size_t i, size_t j, const Trajectory &a,
                           const Trajectory &b, distance_t d) const {
    distance_t start, end;
    std::tie(start, end) = intersection_interval(a[i], d, b[j], b[j + 1]);
    return {start + j, end + j};
  }

  template <typename Trajectory>
  interval get_reachable_b(size_t i, size_t j, const Trajectory &a,
                           const Trajectory &b, distance_t d) const {
    return get_reachable_a(j, i, b, a, d);
  }

  void merge(std::vector<interval> &v, interval i) const {
    if (is_empty_interval(i))
      return;
    if (v.size() && i.first - eps <= v.back().second)
      v.back().second = i.second;
    else
      v.push_back(i);
  }

  template <typename Trajectory>
  distance_t get_last_reachable_point_from_start(const Trajectory &a,
                                                 const Trajectory &b,
                                                 const distance_t d) const {
    size_t j = 0;
    while (j < b.size() - 2 && _dist2(a[0], b[j + 1]) <= sqr(d))
      ++j;

    distance_t result;
    tie(std::ignore, result) = get_reachable_a(0, j, a, b, d);
    return result;
  }

  template <typename Trajectory>
  void get_reachable_intervals(size_t i_min, size_t i_max, size_t j_min,
                               size_t j_max, const curve<Trajectory> &a,
                               const curve<Trajectory> &b, distance_t d,
                               std::vector<interval> &rb,
                               std::vector<interval> &ra,
                               std::vector<interval> &rb_out,
                               std::vector<interval> &ra_out) const {
    interval tb = empty_interval;
    auto it = std::upper_bound(
        rb.begin(), rb.end(),
        interval{j_max, std::numeric_limits<distance_t>::lowest()});
    if (it != rb.begin()) {
      --it;
      if (it->first <= j_max && it->second >= j_min) {
        tb = *it;
      }
    }

    interval ta = empty_interval;
    it = std::upper_bound(
        ra.begin(), ra.end(),
        interval{i_max, std::numeric_limits<distance_t>::lowest()});
    if (it != ra.begin()) {
      --it;
      if (it->first <= i_max && it->second >= i_min) {
        ta = *it;
      }
    }

    if (is_empty_interval(tb) && is_empty_interval(ta))
      return;

    const Trajectory &t1 = a.trajectory();
    const Trajectory &t2 = b.trajectory();
    if (tb.first <= j_min + eps && tb.second >= j_max - eps &&
        ta.first <= i_min + eps && ta.second >= i_max - eps) {
      size_t i_mid = (i_min + 1 + i_max) / 2;
      size_t j_mid = (j_min + 1 + j_max) / 2;
      if (std::sqrt(_dist2(t1[i_mid], t2[j_mid])) +
              std::max(a.curve_length(i_min + 1, i_mid),
                       a.curve_length(i_mid, i_max)) +
              std::max(b.curve_length(j_min + 1, j_mid),
                       b.curve_length(j_mid, j_max)) <=
          d) {
        merge(rb_out, {j_min, j_max});
        merge(ra_out, {i_min, i_max});
        return;
      }
    }

    if (i_min == i_max - 1 && j_min == j_max - 1) {
      interval aa = get_reachable_a(i_max, j_min, t1, t2, d);
      interval bb = get_reachable_b(i_min, j_max, t1, t2, d);

      if (is_empty_interval(ta)) {
        aa.first = std::max(aa.first, tb.first);
      } else if (is_empty_interval(tb)) {
        bb.first = std::max(bb.first, ta.first);
      }

      merge(rb_out, aa);
      merge(ra_out, bb);

    } else {
      if (j_max - j_min > i_max - i_min) {
        std::vector<interval> ra_middle;
        size_t split_position = (j_max + j_min) / 2;
        get_reachable_intervals(i_min, i_max, j_min, split_position, a, b, d,
                                rb, ra, rb_out, ra_middle);
        get_reachable_intervals(i_min, i_max, split_position, j_max, a, b, d,
                                rb, ra_middle, rb_out, ra_out);
      } else {
        std::vector<interval> rb_middle;
        size_t split_position = (i_max + i_min) / 2;
        get_reachable_intervals(i_min, split_position, j_min, j_max, a, b, d,
                                rb, ra, rb_middle, ra_out);
        get_reachable_intervals(split_position, i_max, j_min, j_max, a, b, d,
                                rb_middle, ra, rb_out, ra_out);
      }
    }
  }

public:
  frechet_distance(const squared_distance &dist2 = squared_distance())
      : _dist2(dist2) {}

  frechet_distance(const frechet_distance &) = default;
  frechet_distance(frechet_distance &&) = default;
  frechet_distance &operator=(const frechet_distance &) = default;
  frechet_distance &operator=(frechet_distance &&) = default;

  template <typename Trajectory>
  bool is_bounded_by(const Trajectory &a, const Trajectory &b,
                     distance_t d) const {
    return is_bounded_by(curve<Trajectory>(a, _dist2),
                         curve<Trajectory>(b, _dist2), d);
  }

  template <typename Trajectory>
  bool is_bounded_by(const curve<Trajectory> &c1, const curve<Trajectory> &c2,
                     distance_t d) const {
    assert(c1.size());
    assert(c2.size());

    const Trajectory &t1 = c1.trajectory();
    const Trajectory &t2 = c2.trajectory();
    if (_dist2(t1[0], t2[0]) > sqr(d) ||
        _dist2(t1[t1.size() - 1], t2[t2.size() - 1]) > sqr(d))
      return false;

    if (t1.size() == 1 && t2.size() == 1)
      return true;
    else if (t1.size() == 1)
      return get_dist_to_point_sqr(t2, t1[0]) <= sqr(d);
    else if (t2.size() == 1)
      return get_dist_to_point_sqr(t1, t2[0]) <= sqr(d);

    std::vector<interval> ra, rb, ra_out, rb_out;
    ra.emplace_back(0, get_last_reachable_point_from_start(t1, t2, d));
    rb.emplace_back(0, get_last_reachable_point_from_start(t2, t1, d));

    get_reachable_intervals(0, c1.size() - 1, 0, c2.size() - 1, c1, c2, d, ra,
                            rb, ra_out, rb_out);

    return ra_out.size() &&
           (ra_out.back().second >= c2.size() - static_cast<distance_t>(1.5));
  }
};

/*
* class spatial_index
*/
#include "baldus_bringmann_spatial_index.hpp"

} // bringmanbaldus

} // detail

} // frechetrange

#endif
