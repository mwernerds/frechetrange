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
*         detail            contains implementation details of the three
*                           approaches well-isolated
*             duetschvarhenhold
*                 class FrechetDistance
*                 class Grid
*             submission2
*             submission3
*
*
*/

#include <algorithm> // for std::sort, std::lower_bound, std::upper_bound, and std::max
#include <array>
#include <cmath>      // for std::floor
#include <cstddef>    // for std::size_t
#include <functional> // for std::hash
#include <stdexcept>  // for std::invalid_argument
#include <unordered_map>
#include <utility> // for std::move
#include <vector>

#define ALLOW_FUNCTIONAL

//#include<functional>
#include<cassert>
#include<limits>
#include<iostream>


namespace frechetrange {
namespace detail {
namespace duetschvahrenhold {
// source code from Fabian Dütsch, templatized and generalized by Martin
/**
* Arbitrary constant > 1.0 to signal that a segment is not reachable.
*/
static constexpr double BEGIN_NOT_REACHABLE = 2.0;

template <typename squareddistancefunctional, typename xgetterfunctional,
          typename ygetterfunctional>
class FrechetDistance {
public:
  FrechetDistance(squareddistancefunctional squaredDistance,
                  xgetterfunctional xGetter, ygetterfunctional yGetter)
      : _squaredDistance(squaredDistance), _getX(xGetter), _getY(yGetter),
        _leftSegmentBegins(nullptr), _capacity(0)

  {}
  ~FrechetDistance() { delete[] _leftSegmentBegins; }

  /**
  * Returns whether the Fréchet distance of the passed trajectories is bounded
  * by the passed upper bound.
  * @pre Neither of the trajectories is empty.
  */
  template <typename Trajectory>
  bool isBoundedBy(Trajectory &traj1, Trajectory &traj2,
                    double distanceBound) {
    const double boundSquared = distanceBound * distanceBound;

    // ensure that the corners of the free space diagram are reachable
    if (_squaredDistance(traj1[0], traj2[0]) > boundSquared ||
        _squaredDistance(traj1[traj1.size() - 1], traj2[traj2.size() - 1]) >
            boundSquared)
      return false;

    bool firstIsSmaller = traj1.size() <= traj2.size();
    Trajectory &smallerTraj = firstIsSmaller ? traj1 : traj2;
    Trajectory &biggerTraj = firstIsSmaller ? traj2 : traj1;

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
  squareddistancefunctional _squaredDistance;
  xgetterfunctional _getX;
  ygetterfunctional _getY;

  /**
  * Beginnings of reachable parts of the free space segments on
  * the "frontline", i. e., the right segments of the free space cells lastly
  * processed (and accordingly the left segments of the current column).
  * A segment is reachable, if its respective beginning is <= 1.0.
  */
  double *_leftSegmentBegins;
  /**
  * Allocated size of the array _leftSegmentBegins
  */
  size_t _capacity;

  /**
  * Decides the Fréchet distance problem for a trajectory consisting of only
  * one point.
  */
  template <typename Trajectory, typename point_type>
  bool comparePointToTrajectory(const point_type &p,
                                const Trajectory &trajectory,
                                const double boundSquared) {
    // the first point has already been tested
    for (size_t i = 1; i < trajectory.size(); ++i) {
      if (_squaredDistance(p, trajectory[i]) > boundSquared) {
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
  bool matchInnerPointsMonotonously( Trajectory &points,
                                     Trajectory &segments,
                                    double boundSquared)  {
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
  bool traverseFreeSpaceDiagram( Trajectory &p1, Trajectory &p2,
                                const double boundSquared) {
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
          _squaredDistance(p1[0], p2[colIdx]) > boundSquared) {
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
              _squaredDistance(p1[rowIdx + 1], p2[colIdx + 1]) <= boundSquared;

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
        if (_squaredDistance(p1[lastRow + 1], p2[colIdx + 1]) <= boundSquared &&
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
        if (_squaredDistance(p1[rowIdx + 1], p2[lastColumn + 1]) <=
                boundSquared &&
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
  void clearFrontline(size_t rows) {
    // reserve memory, if necessary
    if (_capacity < rows) {
      delete[] _leftSegmentBegins;
      _leftSegmentBegins = new double[rows];
      _capacity = rows;
    }

    // mark the frontline as not reachable
    std::fill_n(_leftSegmentBegins, rows, BEGIN_NOT_REACHABLE);
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
  template <typename point_type>
  void getLineCircleIntersections(point_type &p1, point_type &p2,
                                   point_type &cp,
                                   double radiusSquared, double &begin,
                                  double &end) {
    // TODO: Adapt to template parameter squareddistancefunctional;
    //       The following assumes the Euclidean distance.

    // Compute the points p1+x*(p2-p1) on the segment from p1 to p2
    // that intersect the circle around cp with radius r:
    // d(cp, p1+x*(p2-p1))^2 = r^2  <=>  a*x^2 + bx + c = 0,
    // where a, b, c, and the auxiliary vectors u and v are defined as follows:

    // u := p2 - p1
    double ux = _getX(p2) - _getX(p1);
    double uy = _getY(p2) - _getY(p1);
    // v := p1 - c
    double vx = _getX(p1) - _getX(cp);
    double vy = _getY(p1) - _getY(cp);

    // a := u^2
    double a = ux * ux + uy * uy;
    // b := 2 * dotProduct(u, v)
    double b = 2 * (vx * ux + vy * uy);
    // c := v^2 - r^2
    double c = vx * vx + vy * vy - radiusSquared;

    // Solve by quadratic formula:
    // x_1,2 = (-b +- sqrt(b^2-4ac)) / 2a
    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
      // The circle does not intersect the line through p1 and p2.
      begin = BEGIN_NOT_REACHABLE;
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
                           double prevOrthogonalSegBegin)  {
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
  *              describedas a scalar value.
  */
  bool isSegmentReachable(double begin) const { return begin <= 1.0; }
};


/*
class Grid
*/
#include "duetsch_vahrenhold_grid.hpp"
} // duetschvahrenhold


namespace bringmanbaldus{








template<typename curvetype, typename pointtype, typename distancetype,  // types
	typename xgetterfunctional, typename ygetterfunctional, typename squareddistancefunctional>     // elementary functions
class FrechetDistance
{

    private:
    static constexpr distancetype eps=10e-5;
    
    typedef curvetype curve;
    typedef distancetype distance_t;
    typedef pointtype point;
//    typedef typename curve::value_type point; @deprecated as we don't want to fix the ::value_type concept for trajectories in details namespace

    // an actual point proxy for having operators
    const xgetterfunctional getx;
    const ygetterfunctional gety;
    squareddistancefunctional dist2;

// arithmetic was in operators, but this was to tedious for an inner class as getx and gety are not meant to be static.
// Therefore, transforming to some sort of assembler.
/*
   DEPRECATED: As we don't know the class point and we don't want to assume anything about it, we should not operate on it. 
   Instead, we have to express everything withing getx and gety and distance_type arithmetics. Not beautiful, but flexible.
   point op_sub(point &a,  point &b)
   {

      getx(a) -= getx(b);
      gety(a) -= gety(b);
      return a;
   }

    point op_minus( point &a,  point &b)
    {
    point tmp(a);
    getx(tmp)  = getx(a)- getx(b);
    gety(tmp) = gety(a) - gety(b);
    return tmp;
    }

   point op_add(point &a,  point &b)
{
    getx(a) += getx(b);
    gety(a) += gety(b);
    return a;
}

    point op_plus( point &a,  point &b)
{
    auto tmp = a;
    op_add(tmp,b);
    return tmp;
}

point op_sdiv(point &a, const distance_t d)
{
    getx(a) /= d;
    gety(a) /= d;
    return a;
}
*/
  template<typename valtype>
    valtype sqr(valtype a) {return a*a;};



distance_t dist_sqr(point& a, point& b)
{
    return(dist2(a,b));
}



    typedef std::pair<distance_t, distance_t> interval; // .first is the startpoint, .second the endpoint (start inclusive, end exclusive)


    



    /*Section 1: special types and their elementary operations*/

     /*point type fun*/
typedef point vec;


    static constexpr interval empty_interval{std::numeric_limits<distance_t>::max(), std::numeric_limits<distance_t>::lowest()};

    /*
    Actually, this would have had to be a free function. But this would pollute the namespace with an interval. we should remove it when done

    std::ostream& operator<<(std::ostream& out, const interval& i)
    {
	out << "[" << i.first << ", " << i.second << "]";
	return out;
    }*/

    inline bool is_empty_interval(const interval& i)
    {
	return i.first >= i.second;
    }

/*    void pout(point & p)
    {
	std::cout << "POINT("<<p[0]<<";" << p[1] << ")" <<std::endl;
	std::cout << "POINT*("<<getx(p)<<";" << gety(p) << ")" <<std::endl;
    }*/
    
interval intersection_interval(point circle_center, distance_t radius, point line_start, point line_end)
{
    // move the circle center to (0, 0) to simplify the calculation
//    line_start -= circle_center;
    distance_t linestartx = getx(line_start)-getx(circle_center);
    distance_t linestarty = gety(line_start)-gety(circle_center);
    distance_t lineendx = getx(line_end)-getx(circle_center);
    distance_t lineendy = gety(line_end)-gety(circle_center);
   

//    op_sub(line_start,circle_center);
//    line_end -= circle_center;
//    op_sub(line_end,circle_center);

    // The line can be represented as line_start + lambda * v

    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 = radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2) + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) = radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 * line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, b := 2 * line_start.x * v.x + 2 * line_start.y * v.y, c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * b + c) = 0
    // <=> lambda^2 + (b / a) * lambda + c / a) = 0
    // <=> lambda1/2 = - (b / 2a) +/- sqrt((b / 2a)^2 - c / a)

     distance_t a = sqr((getx(line_end)-getx(line_start))) + sqr((gety(line_end)-gety(line_start)));
     distance_t b = (linestartx * (lineendx-linestartx) + linestarty * (lineendy-linestarty));
     distance_t c = sqr(linestartx) + sqr(linestarty) - sqr(radius);

     distance_t discriminant = sqr(b / a) - c / a;

    if (discriminant < 0) {
        return empty_interval; // no intersection;
    }

     distance_t lambda1 = - b / a - sqrt(discriminant);
     distance_t lambda2 = - b / a + sqrt(discriminant);

    if (lambda2 < 0 || lambda1 > 1) return empty_interval;
    else return {std::max<distance_t>(lambda1, 0), std::min<distance_t>(lambda2, 1)};
}

    



    /*Section 2: Elementary Frechet operations*/
    template<typename localcurvetype>
distance_t get_dist_to_point_sqr(const localcurvetype& a, point b)
{
    distance_t result = 0;
//    for (point p: a) result = max(result, dist_sqr(p, b));
    return result;
}


inline interval get_reachable_a(size_t i, size_t j, const curve& a, const curve& b, distance_t d)
{
    distance_t start, end;
    std::tie(start, end) = intersection_interval(a[i], d, b[j], b[j + 1]);
    return {start + j, end + j};
}


inline interval get_reachable_b(size_t i, size_t j, const curve& a, const curve& b, distance_t d)
{
    return get_reachable_a(j, i, b, a, d);
}
void merge (std::vector<interval>& v, interval i)
{
    if (is_empty_interval(i)) return;
    if (v.size() && i.first - eps <= v.back().second) v.back().second = i.second;
    else v.push_back(i);
}



distance_t get_last_reachable_point_from_start( curve& a,  curve& b, const distance_t d)
{
    size_t j = 0;
    while (j < b.size() - 2 && dist_sqr(a.front(), b[j + 1]) <= sqr(d)) ++j;
    distance_t result;
    tie(std::ignore, result) = get_reachable_a(0, j, a, b, d);
    return result;
}

// @REMARK: Curvelength has been cached in original in a class curve using a prefix_length generated for each trajectory.
// Omitted now, for easy implementation.

double curve_length(curve&a, size_t i, size_t j)
{
   double d = 0;
   for (auto k = i+1; k <= j; k++)
      d += sqrt(dist2(a[k-1],a[k])); 
   return d;
}


void get_reachable_intervals(size_t i_min, size_t i_max, size_t j_min, size_t j_max, curve& a, curve& b, distance_t d, std::vector<interval>& rb, std::vector<interval>& ra, std::vector<interval>& rb_out, std::vector<interval>& ra_out)
{
	interval tb = empty_interval;
    auto it = std::upper_bound(rb.begin(), rb.end(), interval{j_max, std::numeric_limits<distance_t>::lowest()});
    if (it != rb.begin()) {
        --it;
        if (it->first <= j_max && it->second >= j_min) {
            tb = *it;
        }
    }

    interval ta = empty_interval;
    it = std::upper_bound(ra.begin(), ra.end(), interval{i_max, std::numeric_limits<distance_t>::lowest()});
    if (it != ra.begin()) {
        --it;
        if (it->first <= i_max && it->second >= i_min) {
            ta = *it;
        }
    }

    if (is_empty_interval(tb) && is_empty_interval(ta)) return;

    if (tb.first <= j_min + eps && tb.second >= j_max - eps && ta.first <= i_min + eps && ta.second >= i_max - eps) {
        size_t i_mid = (i_min + 1 + i_max)/2;
		size_t j_mid = (j_min + 1 + j_max)/2;
		if (sqrt(dist2(a[i_mid], b[j_mid])) + std::max(curve_length(a,i_min+1, i_mid),curve_length(a,i_mid, i_max)) + std::max(curve_length(b,j_min+1, j_mid),curve_length(b,j_mid, j_max)) <= d) {
            merge(rb_out, {j_min, j_max});
            merge(ra_out, {i_min, i_max});
            return;
        }
    }

    if (i_min == i_max - 1 && j_min == j_max - 1) {
        interval aa = get_reachable_a(i_max, j_min, a, b, d);
        interval bb = get_reachable_b(i_min, j_max, a, b, d);

        if (is_empty_interval(ta)) {
            aa.first = std::max(aa.first, tb.first);
        }
		else if (is_empty_interval(tb)) { bb.first = std::max(bb.first, ta.first); }

        merge(rb_out, aa);
        merge(ra_out, bb);

    } else {
        if (j_max - j_min > i_max - i_min) {
        	std::vector<interval> ra_middle;
        	size_t split_position = (j_max + j_min) / 2;
        	get_reachable_intervals(i_min, i_max, j_min, split_position, a, b, d, rb, ra, rb_out, ra_middle);
        	get_reachable_intervals(i_min, i_max, split_position, j_max, a, b, d, rb, ra_middle, rb_out, ra_out);
        } else {
        	std::vector<interval> rb_middle;
        	size_t split_position = (i_max + i_min) / 2;
        	get_reachable_intervals(i_min, split_position, j_min, j_max, a, b, d, rb, ra, rb_middle, ra_out);
        	get_reachable_intervals(split_position, i_max, j_min, j_max, a, b, d, rb_middle, ra, rb_out, ra_out);
		}
    }
}


    
   public:

    FrechetDistance (    xgetterfunctional _getx, ygetterfunctional _gety,squareddistancefunctional _dist)
			:getx(_getx),gety(_gety),dist2(_dist){};
   
bool is_frechet_distance_at_most(curve& a, curve& b, distance_t d)
{
    assert(a.size());
    assert(b.size());
    if (dist2(a.front(), b.front()) > d*d || dist2(a.back(), b.back()) > d*d) return false;
    if (a.size() == 1 && b.size() == 1) return true;
    else if (a.size() == 1) return get_dist_to_point_sqr(b, a[0]) <= sqr(d);
    else if (b.size() == 1) return get_dist_to_point_sqr(a, b[0]) <= sqr(d);

    std::vector<interval> ra, rb, ra_out, rb_out;
    ra.push_back({0, get_last_reachable_point_from_start(a, b, d)});
    rb.push_back({0, get_last_reachable_point_from_start(b, a, d)});

    get_reachable_intervals(0, a.size() - 1, 0, b.size() - 1, a, b, d, ra, rb, ra_out, rb_out);

    return ra_out.size() && (ra_out.back().second >= b.size() - static_cast<distance_t>(1.5));

}



   

};






} // bringmanbaldus


} // detail

} // frechetrange

#endif
