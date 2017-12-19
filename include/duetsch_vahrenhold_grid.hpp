#ifndef DUETSCH_VAHRENHOLD_GRID_INC
#define DUETSCH_VAHRENHOLD_GRID_INC


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

  long long toCellNr(double pointCoord) const {
    return static_cast<long long>(std::floor(pointCoord / _meshSize));
  }
  double toCellCoord(double pointCoord) const {
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
  void query(const Trajectory &query, double threshold, Output &output) const {
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

#endif
