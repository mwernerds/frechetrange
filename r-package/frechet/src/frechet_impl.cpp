
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "../../../include/frechetrange.hpp"
#include <array>
#include <functional>
using namespace Rcpp;

struct ContainerAdapterRow {
  const NumericMatrix &m;
  size_t row;
  ContainerAdapterRow(const NumericMatrix &_m, size_t _row)
      : m(_m), row(_row){};
  double operator[](size_t col) const { return m(row, col); }
};

struct ContainerAdapter {

  const NumericMatrix &m;
  typedef ContainerAdapterRow value_type;
  ContainerAdapter(NumericMatrix &_m) : m(_m){};
  ContainerAdapter(const NumericMatrix &_m) : m(_m){};
  const ContainerAdapterRow operator[](size_t row) const {

    return ContainerAdapterRow(m, row);
    ;
  }

  size_t size() const { return m.nrow(); }
};

//@todo: check for 2D in the R interface

auto dist_sqr = [](ContainerAdapterRow r1, ContainerAdapterRow r2) {
  return ((r2[0] - r1[0]) * (r2[0] - r1[0]) +
          (r2[1] - r1[1]) * (r2[1] - r1[1]));
};
auto getx = [](ContainerAdapterRow r1) { return r1[0]; };
auto gety = [](ContainerAdapterRow r1) { return r1[1]; };

frechetrange::detail::duetschvahrenhold::FrechetDistance<
    decltype(dist_sqr), decltype(getx), decltype(gety)>
    fd(dist_sqr, getx, gety);

frechetrange::detail::baldusbringmann::FrechetDistance<
    decltype(dist_sqr), decltype(getx), decltype(gety)>
    fd2(dist_sqr, getx, gety);

// [[Rcpp::export]]
bool internal_frechet_decide_dv(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  ContainerAdapter c1(t1), c2(t2);
  return fd.isBoundedBy(c1, c2, dist);
}

// [[Rcpp::export]]
bool internal_frechet_decide_bb(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  ContainerAdapter c1(t1), c2(t2);
  return fd2.is_frechet_distance_at_most(c1, c2, dist);
}

// --- point types, trajectory types, and global functions ---

constexpr size_t g_DIMENSIONS = 2;
template <size_t dims>
using _point_t = std::array<double, dims>;

auto g_getX = [](const _point_t<g_DIMENSIONS> &p) {
    return p[0];
  };
auto g_getY = [](const _point_t<g_DIMENSIONS> &p) {
    return p[1];
  };
auto g_dist2 = [](const _point_t<g_DIMENSIONS> &p, const _point_t<g_DIMENSIONS> &q) {
    return (q[0] - p[0]) * (q[0] - p[0]) +
           (q[1] - p[1]) * (q[1] - p[1]);
  };

template <size_t dims>
using _trajectory_t = std::vector<_point_t<dims>>;

template <size_t dims>
void _copyMatrixToTrajectory(const NumericMatrix &m, _trajectory_t<dims>& t) {
  for (size_t i = 0; i < m.nrow(); i++) {
    for (size_t d = 0; d < dims; d++)
      t[i][d] = m(i, d);
  }
}

template <size_t dims>
NumericMatrix _trajectoryToMatrix(const _trajectory_t<dims> &t) {
  NumericMatrix m(t.size(), dims);
  for (size_t i = 0; i < t.size(); i++)
    for (size_t d = 0; d < dims; d++)
      m(i, d) = t[i][d];
  return m;
}

template <size_t dims>
List _resultToList(std::vector<const _trajectory_t<dims>*> results) {
  List list(results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    list[i] = _trajectoryToMatrix(*(results[i]));
  }
  return list;
}

// ------------------- Grid -------------------

/// GridDataset caches the data outside of R to have valid container
/// as some of the implementations are not yet compatible.
template <size_t dims> class GridDataset {
public:
  const size_t DIMENSIONS = dims;
  
  typedef frechetrange::detail::duetschvahrenhold::Grid<
      _trajectory_t<dims>, decltype(g_dist2), decltype(g_getX), decltype(g_getY)>
      grid_type;

  GridDataset() : _pGrid(nullptr), _trajectories() {}
  ~GridDataset() { delete _pGrid; }

  grid_type &grid() { return *_pGrid; }
  
  size_t addTrajectory(const NumericMatrix &m) {
    _trajectories.emplace_back(m.nrow());
    _copyMatrixToTrajectory(m, _trajectories.back());
    return (_trajectories.size() - 1);
  }

  size_t size() { return _trajectories.size(); }

  void clear() { _trajectories.clear(); }

  void buildIndex(double meshSize) {
    _pGrid = new grid_type(meshSize, g_dist2, g_getX, g_getY);
    _pGrid->reserve(size());
    for (_trajectory_t<dims> &t : _trajectories)
      _pGrid->insert(std::move(t));
    clear();
    _pGrid->optimize();
  }
  
private:
  grid_type *_pGrid;
  std::vector<_trajectory_t<dims>> _trajectories;
};

std::vector<GridDataset<g_DIMENSIONS>> g_grids;

#define ASSERT_VALID_DATASET(dataset, k)                                       \
  {                                                                            \
    if (k < 0 || k >= dataset.size())                                          \
      throw(std::runtime_error(                                                \
          "Invalid Handle. Create one with createGridDataset"));               \
  }

// [[Rcpp::export]]
size_t internal_createGrid() {
  g_grids.emplace_back();
  return g_grids.size() - 1;
}

// [[Rcpp::export]]
size_t internal_addTrajectoryToGrid(size_t gds, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_grids, gds);
  return g_grids[gds].addTrajectory(m);
}

// [[Rcpp::export]]
size_t internal_clearGrid(size_t gds) {
  ASSERT_VALID_DATASET(g_grids, gds);
  g_grids[gds].clear();
}

// [[Rcpp::export]]
bool internal_buildIndex(size_t gds, double meshSize) {
  ASSERT_VALID_DATASET(g_grids, gds);
  Rcout << "Creating Index from " << g_grids[gds].size() << " trajectories";
  g_grids[gds].buildIndex(meshSize);
}

// [[Rcpp::export]]
List internal_gridRangeQuery(size_t gds, const NumericMatrix &m, double dist,
                             bool materialize = true) {
  ASSERT_VALID_DATASET(g_grids, gds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  
  auto results = g_grids[gds].grid().rangeQuery(t, dist);
  return _resultToList<g_DIMENSIONS>(results);
}

// ------------------- spatial_index -------------------

template <size_t dims>
using _tree_index_type = frechetrange::detail::baldusbringmann::spatial_index<
    _trajectory_t<dims>, decltype(g_dist2), decltype(g_getX), decltype(g_getY)>;

std::vector<_tree_index_type<g_DIMENSIONS>> g_treeIndices;

// [[Rcpp::export]]
size_t internal_createTreeIndex() {
  g_treeIndices.emplace_back(g_dist2, g_getX, g_getY);
  return g_treeIndices.size() - 1;
}

// [[Rcpp::export]]
size_t internal_addTrajectoryToTree(size_t tds, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_treeIndices, tds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  g_treeIndices[tds].add_curve(t);
  return (g_treeIndices[tds].size() - 1);
}

// [[Rcpp::export]]
List internal_treeRangeQuery(size_t tds, const NumericMatrix &m, double dist,
                             bool materialize = true) {
  ASSERT_VALID_DATASET(g_treeIndices, tds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  
  auto results = g_treeIndices[tds].get_close_curves(t, dist);
  return _resultToList<g_DIMENSIONS>(results);
}
