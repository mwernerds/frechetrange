#ifndef DV_FRECHET_DISTANCE_HPP
#define DV_FRECHET_DISTANCE_HPP

#include <algorithm>  // for std::fill_n
#include <cmath>      // for std::sqrt
#include <vector>

#include "../../distance_sqr.hpp"

namespace frechetrange {
namespace detail {
namespace dv {

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
    template <typename trajectory>
    bool is_bounded_by(const trajectory &traj1, const trajectory &traj2,
                       double distance_bound) const {
        const double boundSquared = distance_bound * distance_bound;

        // ensure that the corners of the free space diagram are reachable
        if (_dist2(traj1[0], traj2[0]) > boundSquared ||
            _dist2(traj1[traj1.size() - 1], traj2[traj2.size() - 1]) >
                boundSquared)
            return false;

        bool firstIsSmaller = traj1.size() <= traj2.size();
        const trajectory &smallerTraj = firstIsSmaller ? traj1 : traj2;
        const trajectory &biggerTraj = firstIsSmaller ? traj2 : traj1;

        if (smallerTraj.size() == 1) {
            return comparePointToTrajectory(smallerTraj[0], biggerTraj,
                                            boundSquared);
#ifdef USE_POSITIVE_FILTER
        } else if (positiveFilter(traj1, traj2, boundSquared)) {
            return true;
#endif
        } else if (!matchInnerPointsMonotonously(traj1, traj2, boundSquared) ||
                   !matchInnerPointsMonotonously(traj2, traj1, boundSquared)) {
            // There exists no monotone matchting from one trajectory to the
            // other,
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
    template <typename trajectory, typename point_type>
    bool comparePointToTrajectory(const point_type &p,
                                  const trajectory &t,
                                  const double boundSquared) const {
        // the first point has already been tested
        for (size_t i = 1; i < t.size(); ++i) {
            if (_dist2(p, t[i]) > boundSquared) {
                return false;
            }
        }
        return true;
    }

#ifdef USE_POSITIVE_FILTER
    template <typename trajectory>
    bool positiveFilter(const trajectory &traj1,
                        const trajectory &traj2,
                        const double boundSquared) const {
        // Positive greedy filter as developed by Baldus and Bringmann in
        // "A fast implementation of near neighbors queries for Fréchet distance",
        // SIGSPATIAL'17
        size_t idx1 = 0, idx2 = 0;
        while (idx1 < traj1.size() - 1 && idx2 < traj2.size() - 1) {
            // Distances of three next pairings
            double dist1 = _dist2(traj1[idx1 + 1], traj2[idx2]);
            double dist2 = _dist2(traj1[idx1], traj2[idx2 + 1]);
            double dist12 = _dist2(traj1[idx1 + 1], traj2[idx2 + 1]);

            // Find the minimal distance
            if (dist12 < dist1 && dist12 < dist2) {
                if (dist12 > boundSquared) {
                    return false;
                }

                ++idx1;
                ++idx2;
            } else if (dist1 < dist2) {
                if (dist12 > boundSquared) {
                    return false;
                }

                ++idx1;
            } else {
                if (dist12 > boundSquared) {
                    return false;
                }

                ++idx2;
            }
        }

        // Advance to the end of the first trajectory, if necessary
        while (idx1 < traj1.size() - 2) {
            ++idx1;
            if (_dist2(traj1[idx1], traj2[idx2]) > boundSquared) {
                return false;
            }
        }

        // Advance to the end of the second trajectory, if necessary
        while (idx2 < traj2.size() - 2) {
            ++idx2;
            if (_dist2(traj1[idx1], traj2[idx2]) > boundSquared) {
                return false;
            }
        }

        return true;
    }
#endif

    /**
    * Ensures that the sequence of passed points can be matched to a sequence of
    * points on the passed segments,
    * such that each matching is within the passed distance bound and the
    * sequence of segment points is monotone.
    */
    template <typename trajectory>
    bool matchInnerPointsMonotonously(const trajectory &points,
                                      const trajectory &segments,
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
                                       points[pointIdx], boundSquared, begin,
                                       end);

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
    template <typename trajectory>
    bool traverseFreeSpaceDiagram(const trajectory &p1, const trajectory &p2,
                                  const double boundSquared) const {
        // init the beginnings of reachable parts of the "frontline",
        // i.e., the right segments of the column lastly processed
        // (and accordingly the left segments of the current column)
        clearFrontline(p1.size() - 1);

        // the bottom-left free space corner has been ensured to be reachable
        _leftSegmentBegins[0] = 0.0;
        // beginning of the reachable part of the bottommost segment of the
        // current
        // column described as a scalar value
        double bottommostSegmentBegin = 0.0;
        // bottommost reachable segment of the current column
        size_t bottommostReachableRow = 0;

        const size_t lastRow = p1.size() - 2;
        const size_t lastColumn = p2.size() - 2;
        // compute the reachable parts of the free space diagram's cells
        // columnwise
        for (size_t colIdx = 0; colIdx < lastColumn; ++colIdx) {
            // compute whether the current column's bottommost segment is
            // reachable
            if (bottommostSegmentBegin == 0.0 &&
                _dist2(p1[0], p2[colIdx]) > boundSquared) {
                bottommostSegmentBegin = BEGIN_NOT_REACHABLE;
            }

            // beginning of the reachable part of the current bottom segment
            double currBottomSegBegin = bottommostSegmentBegin;
            double currBottomSegEnd = currBottomSegBegin;

            double currRightSegBegin, currRightSegEnd;
            // traverse this columns' reachable cells
            for (size_t rowIdx = bottommostReachableRow; rowIdx < lastRow;
                 ++rowIdx) {
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
                        getLineCircleIntersections(
                            p1[rowIdx], p1[rowIdx + 1], p2[colIdx + 1],
                            boundSquared, currRightSegBegin, currRightSegEnd);
                        currRightSegBegin = getReachableBegin(
                            currRightSegBegin, currRightSegEnd,
                            _leftSegmentBegins[rowIdx], currBottomSegBegin);
                    }

                    // compute the reachable part of the current cell's top
                    // segment
                    double currTopSegBegin, currTopSegEnd;
                    if (isTopRightFree &&
                        _leftSegmentBegins[rowIdx + 1] <= 0.0) {
                        // the entire top segment is reachable
                        currTopSegBegin = 0.0;
                        currTopSegEnd = 1.0;
                    } else {
                        getLineCircleIntersections(
                            p2[colIdx], p2[colIdx + 1], p1[rowIdx + 1],
                            boundSquared, currTopSegBegin, currTopSegEnd);
                        currTopSegBegin = getReachableBegin(
                            currTopSegBegin, currTopSegEnd, currBottomSegBegin,
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
                    getLineCircleIntersections(
                        p1[lastRow], p1[lastRow + 1], p2[colIdx + 1],
                        boundSquared, currRightSegBegin, currRightSegEnd);
                    currRightSegBegin = getReachableBegin(
                        currRightSegBegin, currRightSegEnd,
                        _leftSegmentBegins[lastRow], currBottomSegBegin);
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
        for (size_t rowIdx = bottommostReachableRow; rowIdx < lastRow;
             ++rowIdx) {
            // ensure that this cell is reachable
            if (isSegmentReachable(_leftSegmentBegins[rowIdx]) ||
                isSegmentReachable(currBottomSegBegin)) {
                // compute the reachable part of the current cell's top segment
                double currTopSegBegin, currTopSegEnd;
                if (_dist2(p1[rowIdx + 1], p2[lastColumn + 1]) <=
                        boundSquared &&
                    _leftSegmentBegins[rowIdx + 1] <= 0.0) {
                    // the entire top segment is reachable
                    currBottomSegBegin = 0.0;
                } else {
                    getLineCircleIntersections(
                        p2[lastColumn], p2[lastColumn + 1], p1[rowIdx + 1],
                        boundSquared, currTopSegBegin, currTopSegEnd);
                    currBottomSegBegin = getReachableBegin(
                        currTopSegBegin, currTopSegEnd, currBottomSegBegin,
                        _leftSegmentBegins[rowIdx]);
                }
            }
        }

        // return whether the top-right corner is reachable via its bottom
        // segment
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
        // where a, b, c, and the auxiliary vectors u and v are defined as
        // follows:
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
            return (currSegmentEnd >= prevParallelSegBegin)
                       ? prevParallelSegBegin
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

}  // namespace dv
}  // namespace detail
}  // namespace frechetrange

#endif
