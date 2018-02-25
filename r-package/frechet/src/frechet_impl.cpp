#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include <array>
#include <functional>

#include "../../../include/frechetrange.hpp"

using Rcpp::NumericMatrix;
using Rcpp::List;
using Rcpp::Rcout;

typedef NumericMatrix::ConstRow point_adapter;

struct trajectory_adapter {
  const NumericMatrix &_m;
  
  typedef point_adapter value_type;

  trajectory_adapter(const NumericMatrix &m) : _m(m){};
  point_adapter operator[](size_t idx) const {
    return _m.row(idx);
  }

  size_t size() const { return _m.nrow(); }
};

//@todo: check for 2D in the R interface

struct get_adapter_coord {
  template <size_t dim> static double get(const point_adapter &c) {
    return c[dim];
  }
};

frechetrange::detail::duetschvahrenhold::frechet_distance<
    2, get_adapter_coord> fd;

frechetrange::detail::baldusbringmann::frechet_distance<
    2, get_adapter_coord> fd2;

// [[Rcpp::export]]
bool internal_frechet_decide_dv(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  trajectory_adapter c1(t1), c2(t2);
  return fd.is_bounded_by(c1, c2, dist);
}

// [[Rcpp::export]]
bool internal_frechet_decide_bb(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  trajectory_adapter c1(t1), c2(t2);
  return fd2.is_bounded_by(c1, c2, dist);
}

// --- point types, trajectory types, and global functions ---

constexpr size_t g_DIMENSIONS = 2;
template <size_t dims>
using _point_t = std::array<double, dims>;

struct get_point_coord {
  template <size_t dim> static double get(const _point_t<g_DIMENSIONS> &p) {
    return std::get<dim>(p);
  }
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

/// grid_data caches the data outside of R to have valid container
/// as some of the implementations are not yet compatible.
template <size_t dims> class grid_data {
public:
  const size_t DIMENSIONS = dims;
  
  typedef frechetrange::detail::duetschvahrenhold::grid<dims,
      _trajectory_t<dims>, get_point_coord>
      grid_type;

  grid_data() : _pGrid(nullptr), _trajectories() {}
  ~grid_data() { delete _pGrid; }

  grid_type &grid() { return *_pGrid; }
  
  size_t add_trajectory(const NumericMatrix &m) {
    _trajectories.emplace_back(m.nrow());
    _copyMatrixToTrajectory(m, _trajectories.back());
    return (_trajectories.size() - 1);
  }

  size_t size() { return _trajectories.size(); }

  void clear() { _trajectories.clear(); }

  void build_index(double meshSize) {
    _pGrid = new grid_type(meshSize);
    _pGrid->reserve(size());
    for (_trajectory_t<dims> &t : _trajectories)
      _pGrid->insert(std::move(t));
    clear();
    _pGrid->build_index();
  }
  
private:
  grid_type *_pGrid;
  std::vector<_trajectory_t<dims>> _trajectories;
};

std::vector<grid_data<g_DIMENSIONS>> g_grids;

#define ASSERT_VALID_DATASET(dataset, k)                                       \
  {                                                                            \
    if (k < 0 || k >= dataset.size())                                          \
      throw(std::runtime_error(                                                \
          "Invalid Handle."));               \
  }

// [[Rcpp::export]]
size_t internal_createGrid() {
  g_grids.emplace_back();
  return g_grids.size() - 1;
}

// [[Rcpp::export]]
size_t internal_addTrajectoryToGrid(size_t gds, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_grids, gds);
  return g_grids[gds].add_trajectory(m);
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
  g_grids[gds].build_index(meshSize);
}

// [[Rcpp::export]]
List internal_gridRangeQuery(size_t gds, const NumericMatrix &m, double dist,
                             bool materialize = true) {
  ASSERT_VALID_DATASET(g_grids, gds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  
  auto results = g_grids[gds].grid().range_query(t, dist);
  return _resultToList<g_DIMENSIONS>(results);
}

// ------------------- spatial_index -------------------

template <size_t dims>
using _tree_index_type = frechetrange::detail::baldusbringmann::spatial_index<
    2, _trajectory_t<dims>, get_point_coord>;

std::vector<_tree_index_type<g_DIMENSIONS>> g_treeIndices;

// [[Rcpp::export]]
size_t internal_createTreeIndex() {
  g_treeIndices.emplace_back();
  return g_treeIndices.size() - 1;
}

// [[Rcpp::export]]
size_t internal_addTrajectoryToTree(size_t tds, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_treeIndices, tds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  g_treeIndices[tds].insert(t);
  return (g_treeIndices[tds].size() - 1);
}

// [[Rcpp::export]]
List internal_treeRangeQuery(size_t tds, const NumericMatrix &m, double dist,
                             bool materialize = true) {
  ASSERT_VALID_DATASET(g_treeIndices, tds);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  
  auto results = g_treeIndices[tds].range_query(t, dist);
  return _resultToList<g_DIMENSIONS>(results);
}
