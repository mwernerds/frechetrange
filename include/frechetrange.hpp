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
  bool isBoundedBy(const Trajectory &traj1, const Trajectory &traj2,
                   const double distanceBound) {
    const double boundSquared = distanceBound * distanceBound;

    // ensure that the corners of the free space diagram are reachable
    if (_squaredDistance(traj1[0], traj2[0]) > boundSquared ||
        _squaredDistance(traj1[traj1.size() - 1], traj2[traj2.size() - 1]) >
            boundSquared)
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
  void getLineCircleIntersections(const point_type &p1, const point_type &p2,
                                  const point_type &cp,
                                  const double radiusSquared, double &begin,
                                  double &end) const {
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
  *              describedas a scalar value.
  */
  bool isSegmentReachable(double begin) const { return begin <= 1.0; }
};

/**
* A grid of fixed mesh size spanning the Euclidean plane.
* It stores dataset trajectories in its cells to
* efficiently perform range queries on their MBR corners.
* Trajectories are mapped to cells using one of the four bounding box corners.
*/
template <typename Trajectory, typename squareddistancefunctional,
          typename xgetterfunctional, typename ygetterfunctional>
class Grid {
public:
  /**
  * Creates a grid of the specified mesh size.
  */
  Grid(double meshSize, squareddistancefunctional squaredDistance,
       xgetterfunctional xGetter, ygetterfunctional yGetter)
      : _meshSize(meshSize), _maps(), _expectedQueryCost{{0, 0, 0, 0}},
        _useLeftBorder(), _useBottomBorder(), _optimized(false),
        _decider(squaredDistance, xGetter, yGetter), _getX(xGetter),
        _getY(yGetter) {}
  Grid(const Grid &) = default;
  Grid(Grid &&) = default;
  Grid &operator=(const Grid &) = default;
  Grid &operator=(Grid &&) = default;

  /**
  * Returns the mesh size of this grid.
  */
  double getMeshSize() const { return _meshSize; }

  /**
  * Reserves internal storage so that the indicated number
  * of inserts can be performed without resizing.
  */
  void reserve(size_t numTrajectories) {
    for (size_t i = 0; i < NUM_MAPS; ++i) {
      if (!_optimized || i == toMapIndex(_useLeftBorder, _useBottomBorder))
        _maps[i].reserve(numTrajectories);
    }
  }

  void insert(const Trajectory &trajectory) {
    if (trajectory.size() > 0)
      insertImpl(MBR(trajectory, _getX, _getY));
  }

  void insert(Trajectory &&trajectory) {
    if (trajectory.size() > 0)
      insertImpl(MBR(std::move(trajectory), _getX, _getY));
  }

  /**
  * Determines which MBR corner to use during hashing to minimize
  * the expected query cost. Furthermore, the grid cells are sorted
  * to achieve better query times.
  */
  void optimize() {
    if (_optimized)
      return;

    // only keep the grid with the best expected performance
    chooseBestMap();
    for (size_t i = 0; i < NUM_MAPS; ++i) {
      if (i != toMapIndex(_useLeftBorder, _useBottomBorder))
        _maps[i].clear();
    }

    // sort the grid cells
    if (_useLeftBorder)
      if (_useBottomBorder)
        sortCells<true, true>();
      else
        sortCells<true, false>();
    else if (_useBottomBorder)
      sortCells<false, true>();
    else
      sortCells<false, false>();

    _optimized = true;
  }

  std::vector<const Trajectory *> rangeQuery(const Trajectory &query,
                                             double distanceThreshold) {
    std::vector<const Trajectory *> resultSet;
    auto pushBackResult = [&resultSet](const Trajectory *t) {
      resultSet.push_back(t);
    };
    rangeQuery(query, distanceThreshold, pushBackResult);
    return resultSet;
  }

  /**
  * @param output Supports the method operator()(const Trajectory *)
  *               to output the result trajectories.
  */
  template <typename OutputFunctional>
  void rangeQuery(const Trajectory &query, double distanceThreshold,
                  OutputFunctional &output) {
    if (query.size() == 0)
      return;
    else if (distanceThreshold > _meshSize)
      throw std::invalid_argument(
          "The distance threshold is grater than the mesh size.");

    if (!_optimized)
      chooseBestMap();

    if (_useLeftBorder)
      if (_useBottomBorder)
        this->template query<true, true, OutputFunctional>(
            query, distanceThreshold, output);
      else
        this->template query<true, false, OutputFunctional>(
            query, distanceThreshold, output);
    else if (_useBottomBorder)
      this->template query<false, true, OutputFunctional>(
          query, distanceThreshold, output);
    else
      this->template query<false, false, OutputFunctional>(
          query, distanceThreshold, output);
  }

private:
  /**
  * Integer coordinates of a grid cell
  */
  struct CellNr {
    long long x, y;
    bool operator==(const CellNr &other) const {
      return x == other.x && y == other.y;
    }
  };

  struct MBR {
    // TODO: evaluate using pointer instead
    Trajectory trajectory;
    /**
    * Coordinates of the borders of the bounding box
    */
    double minX, maxX, minY, maxY;

    MBR(const Trajectory &t, const xgetterfunctional &getX,
        const ygetterfunctional &getY)
        : trajectory(t) {
      initBorders(getX, getY);
    }
    MBR(Trajectory &&t, const xgetterfunctional getX,
        const ygetterfunctional &getY)
        : trajectory(std::move(t)) {
      initBorders(getX, getY);
    }
    MBR() = default;
    MBR(const MBR &) = default;
    MBR(MBR &&) = default;
    MBR &operator=(const MBR &) = default;
    MBR &operator=(MBR &) = default;

    void initBorders(const xgetterfunctional &getX,
                     const ygetterfunctional &getY) {
      minX = maxX = getX(trajectory[0]);
      minY = maxY = getY(trajectory[0]);

      for (size_t i = 1; i < trajectory.size(); ++i) {
        auto xCoord = getX(trajectory[i]);
        if (xCoord < minX)
          minX = xCoord;
        else if (xCoord > maxX)
          maxX = xCoord;

        auto yCoord = getY(trajectory[i]);
        if (yCoord < minY)
          minY = yCoord;
        else if (yCoord > maxY)
          maxY = yCoord;
      }
    }
    /**
    * Returns the coordinate of the specified border of the bounding box.
    * @tparam xDim Whether a x- or y-coordinate is returned
    * @tparam first Whether the first (i. e., left resp. bottom) or
    *               second (i.e., right resp. top) border is returned
    */
    template <bool xDim, bool first> double getBorder() const {
      // the compiler hopefully eliminates the branches
      if (xDim)
        if (first)
          return minX;
        else
          return maxX;
      else if (first)
        return minY;
      else
        return maxY;
    }
  };

  /**
  * Compare MBRs by the specified border coordinate.
  */
  template <bool xDim, bool first> struct MBRComparator {
    bool operator()(const MBR &m1, const MBR &m2) const {
      return m1.template getBorder<xDim, first>() <
             m2.template getBorder<xDim, first>();
    }
  };

  /**
  * Cell of the grid
  */
  struct Cell {
    using Bucket = std::vector<MBR>;
    /**
    * Whether this cell is sorted by x- or y-coordinates
    */
    bool xSorted;
    /**
    * List of dataset trajectories stored in this grid cell
    */
    Bucket bucket;

    Cell() = default;
    Cell(const Cell &) = default;
    Cell(Cell &&) = default;
    Cell &operator=(const Cell &) = default;
    Cell &operator=(Cell &) = default;

    template <bool left, bool bottom> void sort() {
      // decide whether to sort by x- or y-coordinates
      xSorted = chooseSortingOrder<left, bottom>();
      if (xSorted)
        std::sort(bucket.begin(), bucket.end(), MBRComparator<true, left>());
      else
        std::sort(bucket.begin(), bucket.end(), MBRComparator<false, bottom>());
    }
    template <bool left, bool bottom> bool chooseSortingOrder() {
      // find extremal MBR corner coordinates
      double minX = bucket[0].template getBorder<true, left>();
      double maxX = minX;
      double minY = bucket[0].template getBorder<false, bottom>();
      double maxY = minY;

      for (size_t i = 1; i < bucket.size(); ++i) {
        double x = bucket[i].template getBorder<true, left>();
        if (x < minX) {
          minX = x;
        } else if (x > maxX) {
          maxX = x;
        }
        double y = bucket[i].template getBorder<false, bottom>();
        if (y < minY) {
          minY = y;
        } else if (y > maxY) {
          maxY = y;
        }
      }

      // choose the dimension with greater diffusion as sorting order
      double xRange = maxX - minX;
      double yRange = maxY - minY;
      return xRange >= yRange;
    }
  };

  struct CellHasher {
    size_t operator()(const CellNr &cell) const {
      std::hash<long long> longHasher;
      size_t hash = longHasher(cell.x) + 0x9e3779b9;
      hash ^= longHasher(cell.y) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      return hash;
    }
  };

  /**
  * Minimal number of containing trajectories for a cell to be sorted.
  */
  static constexpr size_t MIN_SORT_SIZE = 16;
  /**
  * Describes which horizontal resp. vertical cell borders are crossed
  * by a range query
  */
  using Crossings = unsigned int;
  /**
  * No cell border is crossed
  */
  static constexpr Crossings CROSSES_NONE = 0;
  /**
  * The left resp. bottom cell border is crossed
  */
  static constexpr Crossings CROSSES_BEGIN = 2;
  /**
  * The right resp. top cell border is crossed
  */
  static constexpr Crossings CROSSES_END = 1;

  /**
  * Mesh size of this grid, i. e., the side length of cells
  */
  double _meshSize;
  using Map = std::unordered_map<CellNr, Cell, CellHasher>;
  static constexpr size_t NUM_MAPS = 4;
  /**
  * Map from the 2D coordinates to cells of dataset trajectories
  * for each of the four MBR corners.
  */
  std::array<Map, NUM_MAPS> _maps;
  /**
  * Sum of squared cell sizes of each mapping.
  */
  std::array<size_t, NUM_MAPS> _expectedQueryCost;
  /**
  * Whether the bounding box corner, according to which queries are mapped to
  * grid cells, lies on the left (or right) border of the bounding box.
  */
  bool _useLeftBorder;
  /**
  * Whether the bounding box corner, according to which queries are mapped to
  * grid cells, lies on the bottom (or top) border of the bounding box.
  */
  bool _useBottomBorder;
  /**
  * Whether the cells have been sorted.
  */
  bool _optimized;

  mutable FrechetDistance<squareddistancefunctional, xgetterfunctional,
                          ygetterfunctional>
      _decider;
  squareddistancefunctional _squaredDistance;
  xgetterfunctional _getX;
  ygetterfunctional _getY;

  long long toCellNr(double pointCoord) {
    return static_cast<long long>(std::floor(pointCoord / _meshSize));
  }
  double toCellCoord(double pointCoord) {
    return std::floor(pointCoord / _meshSize) * _meshSize;
  }

  size_t toMapIndex(bool left, bool bottom) const {
    return 2 * static_cast<size_t>(left) + static_cast<size_t>(bottom);
  }

  void insertImpl(MBR &&mbr) {
    if (!_optimized) {
      insertInMap<false, false>(MBR(mbr));
      insertInMap<false, true>(MBR(mbr));
      insertInMap<true, false>(MBR(mbr));
      insertInMap<true, true>(std::move(mbr));
    } else {
      if (_useLeftBorder)
        if (_useBottomBorder)
          insertInMap<true, true>(std::move(MBR(mbr)));
        else
          insertInMap<true, false>(std::move(MBR(mbr)));
      else if (_useBottomBorder)
        insertInMap<false, true>(std::move(MBR(mbr)));
      else
        insertInMap<false, false>(std::move(MBR(mbr)));
    }
  }

  template <bool left, bool bottom> void insertInMap(MBR &&mbr) {
    // integer coordinates of the grid cell
    CellNr cellNr{toCellNr(mbr.template getBorder<true, left>()),
                  toCellNr(mbr.template getBorder<false, bottom>())};
    Cell &cell = _maps[toMapIndex(left, bottom)][cellNr];
    typename Cell::Bucket &bucket = cell.bucket;

    if (!_optimized) {
      // append the trajectory to the cell's bucket
      bucket.emplace_back(std::move(mbr));
      // update the sum of squared bucket sizes:
      // (b + 1)^2 = b^2 + 2*(b+1) - 1
      _expectedQueryCost[toMapIndex(left, bottom)] += 2 * bucket.size() - 1;
    } else {
      // insert in the sorted bucket
      auto insertPos = cell.xSorted
                           ? std::upper_bound(bucket.begin(), bucket.end(), mbr,
                                              MBRComparator<true, left>())
                           : std::upper_bound(bucket.begin(), bucket.end(), mbr,
                                              MBRComparator<false, bottom>());
      bucket.emplace(insertPos, std::move(mbr));
    }
  }

  void chooseBestMap() {
    // find the map with best expected query cost
    size_t lowestCost = _expectedQueryCost[0];
    _useLeftBorder = _useBottomBorder = false;

    if (_expectedQueryCost[1] < lowestCost) {
      lowestCost = _expectedQueryCost[1];
      _useLeftBorder = false;
      _useBottomBorder = true;
    }

    if (_expectedQueryCost[2] < lowestCost) {
      lowestCost = _expectedQueryCost[2];
      _useLeftBorder = true;
      _useBottomBorder = false;
    }

    if (_expectedQueryCost[3] < lowestCost) {
      lowestCost = _expectedQueryCost[3];
      _useLeftBorder = true;
      _useBottomBorder = true;
    }
  }

  template <bool left, bool bottom> void sortCells() {
    for (auto &iter : _maps[toMapIndex(left, bottom)]) {
      Cell &cell = iter.second;
      if (cell.bucket.size() >= MIN_SORT_SIZE)
        cell.template sort<left, bottom>();
    }
  }

  template <bool left, bool bottom, typename Output>
  void query(const Trajectory &query, double threshold, Output &output) {
    MBR queryMBR(query, _getX, _getY);
    // check which horizontal neighbor cells need to be visited
    double cellCoordX = toCellCoord(queryMBR.template getBorder<true, left>());
    bool visitLeft =
        queryMBR.template getBorder<true, left>() - threshold < cellCoordX;
    bool visitRight = queryMBR.template getBorder<true, left>() + threshold >=
                      cellCoordX + _meshSize;

    // check which vertical neighbor cells need to be visited
    double cellCoordY =
        toCellCoord(queryMBR.template getBorder<false, bottom>());
    bool visitBottom =
        queryMBR.template getBorder<false, bottom>() - threshold < cellCoordY;
    bool visitTop = queryMBR.template getBorder<false, bottom>() + threshold >=
                    cellCoordY + _meshSize;

    // memorize the crossed cell borders
    Crossings crossedVerticals =
        static_cast<Crossings>(visitLeft) * CROSSES_BEGIN +
        static_cast<Crossings>(visitRight) * CROSSES_END;
    Crossings crossedHorizontals =
        static_cast<Crossings>(visitBottom) * CROSSES_BEGIN +
        static_cast<Crossings>(visitTop) * CROSSES_END;

    // integral coordinates of the cell the query trajectory is mapped to
    CellNr cellNr{toCellNr(queryMBR.template getBorder<true, left>()),
                  toCellNr(queryMBR.template getBorder<false, bottom>())};

    // visit the center cell
    checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                    crossedVerticals, crossedHorizontals,
                                    output);
    // visit the bottom cell
    --cellNr.y;
    if (visitBottom)
      checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                      crossedVerticals, CROSSES_END, output);
    // visit the top cell
    cellNr.y += 2;
    if (visitTop)
      checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                      crossedVerticals, CROSSES_BEGIN, output);

    --cellNr.x;
    if (visitLeft) {
      // visit the top-left cell
      if (visitTop)
        checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                        CROSSES_END, CROSSES_BEGIN, output);
      // visit the left cell
      --cellNr.y;
      checkCell<left, bottom, Output>(cellNr, queryMBR, threshold, CROSSES_END,
                                      crossedHorizontals, output);
      // visit the bottom-left cell
      --cellNr.y;
      if (visitBottom)
        checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                        CROSSES_END, CROSSES_END, output);
      cellNr.y += 2;
    }

    cellNr.x += 2;
    if (visitRight) {
      // visit the top-right cell
      if (visitTop)
        checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                        CROSSES_BEGIN, CROSSES_BEGIN, output);
      // visit the right cell
      --cellNr.y;
      checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                      CROSSES_BEGIN, crossedHorizontals,
                                      output);
      // visit the bottom-right cell
      --cellNr.y;
      if (visitBottom)
        checkCell<left, bottom, Output>(cellNr, queryMBR, threshold,
                                        CROSSES_BEGIN, CROSSES_END, output);
    }
  }

  template <bool left, bool bottom, typename Output>
  void checkCell(const CellNr &cellNr, const MBR &queryMBR, double threshold,
                 Crossings crossedVerticals, Crossings crossedHorizontals,
                 Output &output) const {
    const Map &map = _maps[toMapIndex(left, bottom)];
    // ensure that the cell containts some trajectories
    typename Map::const_iterator iter = map.find(cellNr);
    if (iter == map.end())
      return;

    const Cell &cell = iter->second;
    // traverse the contained trajectories according to the sorting order
    if (cell.xSorted) {
      this->template traverseBucket<true, left, Output>(
          cell.bucket, queryMBR, threshold, crossedVerticals, output);
    } else {
      this->template traverseBucket<false, bottom, Output>(
          cell.bucket, queryMBR, threshold, crossedHorizontals, output);
    }
  }

  template <bool xDim, bool firstBorder, typename Output>
  void traverseBucket(const typename Cell::Bucket &bucket, const MBR &queryMBR,
                      double threshold, Crossings crossedBorders,
                      Output &output) const {
    if (bucket.size() < MIN_SORT_SIZE || !_optimized) {
      // check each trajectory of this cell,
      // as they are not sorted
      for (const MBR &trajectory : bucket) {
        checkTrajectory<false, false, false>(queryMBR, threshold, trajectory,
                                             output);
      }
    } else { // the trajectories are sorted
      // choose the traversing order and beginning depending on
      // which cell borders are crossed
      if (crossedBorders == CROSSES_END) {
        // traverse the trajectories backwards,
        // until the beginning of the active range is reached
        double activeRangeBegin =
            queryMBR.template getBorder<xDim, firstBorder>() - threshold;
        for (size_t i = bucket.size() - 1;
             bucket[i].template getBorder<xDim, firstBorder>() >=
             activeRangeBegin;
             --i) {
          checkTrajectory<true, xDim, firstBorder>(queryMBR, threshold,
                                                   bucket[i], output);

          if (i == 0) {
            break;
          }
        }
      } else {
        // traverse the trajectories forwards,
        // until the end of the active range is reached
        auto searchBegin = bucket.cbegin();
        auto searchEnd = bucket.cend();

        // search for the beginning of the active range,
        // if the range query does not span the first,
        // i. e., left or bottom, cell border
        if (crossedBorders == CROSSES_NONE) {
          double activeRangeBegin =
              queryMBR.template getBorder<xDim, firstBorder>() - threshold;
          // decide whether to find the beginning of the active range
          // by binary searching or linear searching
          if (useBinarySearch(bucket.size(), threshold)) {
            // binary search for the beginning of the active range
            searchBegin = std::lower_bound(
                searchBegin, searchEnd, activeRangeBegin,
                [](const MBR &m, double d) {
                  return m.template getBorder<xDim, firstBorder>() < d;
                });
          } else {
            // linear search for the beginning of the active range
            while (searchBegin != searchEnd &&
                   searchBegin->template getBorder<xDim, firstBorder>() <
                       activeRangeBegin) {
              ++searchBegin;
            }
          }
        }

        // traverse the trajectories forwards,
        // until the end of the active range is reached
        double activeRangeEnd =
            queryMBR.template getBorder<xDim, firstBorder>() + threshold;
        while (searchBegin != searchEnd &&
               searchBegin->template getBorder<xDim, firstBorder>() <=
                   activeRangeEnd) {
          checkTrajectory<true, xDim, firstBorder>(queryMBR, threshold,
                                                   *searchBegin, output);
          ++searchBegin;
        }
      }
    }
  }

  template <bool prechecked, bool xDim, bool firstBorder,
            typename OutputFunctional>
  void checkTrajectory(const MBR &queryMBR, double threshold,
                       const MBR &trajMBR, OutputFunctional &output) const {
    // ensure that all bounding box borders are within range
    if (mbrsWithinRange<prechecked, xDim, firstBorder>(queryMBR, threshold,
                                                       trajMBR)) {
      // append the dataset trajectory to the output,
      // if it is within Fréchet distance of the query trajectory
      if (squaredfarthestBBDistance(queryMBR, trajMBR) <=
              threshold * threshold ||
          _decider.template isBoundedBy<Trajectory>(
              queryMBR.trajectory, trajMBR.trajectory, threshold)) {
        output(&(trajMBR.trajectory));
      }
    }
  }

  template <bool prechecked, bool xDim, bool firstBorder>
  bool mbrsWithinRange(const MBR &queryMBR, double threshold,
                       const MBR &trajMBR) const {
    if (prechecked) {
      // ensure that the three bounding box borders that have not been
      // checked, yet, are within range
      return bordersWithinRange<!xDim, !firstBorder>(queryMBR, threshold,
                                                     trajMBR) &&
             bordersWithinRange<xDim, !firstBorder>(queryMBR, threshold,
                                                    trajMBR) &&
             bordersWithinRange<!xDim, !firstBorder>(queryMBR, threshold,
                                                     trajMBR);
    } else {
      // ensure that all bounding box borders are within range
      return bordersWithinRange<true, true>(queryMBR, threshold, trajMBR) &&
             bordersWithinRange<true, false>(queryMBR, threshold, trajMBR) &&
             bordersWithinRange<false, true>(queryMBR, threshold, trajMBR) &&
             bordersWithinRange<false, false>(queryMBR, threshold, trajMBR);
    }
  }

  template <bool xDim, bool firstBorder>
  bool bordersWithinRange(const MBR &queryMBR, double threshold,
                          const MBR &trajMBR) const {
    return queryMBR.template getBorder<xDim, firstBorder>() - threshold <=
               trajMBR.template getBorder<xDim, firstBorder>() &&
           trajMBR.template getBorder<xDim, firstBorder>() <=
               queryMBR.template getBorder<xDim, firstBorder>() + threshold;
  }

  double squaredfarthestBBDistance(const MBR &queryMBR,
                                   const MBR &trajMBR) const {
    double dx1 = trajMBR.template getBorder<true, false>() -
                 queryMBR.template getBorder<true, true>();
    double dx2 = queryMBR.template getBorder<true, false>() -
                 trajMBR.template getBorder<true, true>();
    double dy1 = trajMBR.template getBorder<false, false>() -
                 queryMBR.template getBorder<false, true>();
    double dy2 = queryMBR.template getBorder<false, false>() -
                 trajMBR.template getBorder<false, true>();
    return std::max(dx1 * dx1, dx2 * dx2) + std::max(dy1 * dy1, dy2 * dy2);
  }

  bool useBinarySearch(size_t bucketSize, double threshold) const {
    // minimal number of elements to choose binary over linear search
    constexpr size_t MIN_BINARY_SEARCH_SIZE = 32;
    // expected ratio of elements of this bucket to skip;
    // it is not negative, as _meshSize > 2*threshold
    // otherwise, a cell border would be crossed and
    // this method would not be called
    double expectedRatioElemsToSkip = 0.5 - threshold / _meshSize;
    return expectedRatioElemsToSkip * bucketSize >= MIN_BINARY_SEARCH_SIZE;
  }
};
} // duetschvahrenhold


namespace bringmanbaldus{


/*THIS SECTION IS IN BRANCH mwernerds/bb, as it is unexpectedly difficult to resolve. The branch is not yet online...*/





} // bringmanbaldus


} // detail

} // frechetrange

#endif
