#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include <array>
#include <stdexcept> // for std::runtime_error
#include <utility>   // for std::move
#include <vector>

//#define ENABLE_MULTITHREADING
#include "../../../include/frechetrange/frechetrange.hpp"

using Rcpp::NumericMatrix;
using Rcpp::List;
using Rcpp::Rcout;

typedef NumericMatrix::ConstRow point_adapter;

struct trajectory_adapter {
  const NumericMatrix &_m;

  typedef point_adapter value_type;

  trajectory_adapter(const NumericMatrix &m) : _m(m){};
  point_adapter operator[](size_t idx) const { return _m.row(idx); }

  size_t size() const { return _m.nrow(); }
};

//@todo: check for 2D in the R interface

struct get_adapter_coord {
  template <size_t dim> static double get(const point_adapter &c) {
    return c[dim];
  }
};

frechetrange::detail::dv::frechet_distance<2, get_adapter_coord>
    fd;

frechetrange::detail::bb::frechet_distance<2, get_adapter_coord>
    fd2;

// [[Rcpp::export]]
bool internal_dv_frechet_decide(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  trajectory_adapter c1(t1), c2(t2);
  return fd.is_bounded_by(c1, c2, dist);
}

// [[Rcpp::export]]
bool internal_bb_frechet_decide(NumericMatrix &t1, NumericMatrix &t2,
                                double dist) {
  trajectory_adapter c1(t1), c2(t2);
  return fd2.is_bounded_by(c1, c2, dist);
}

// --- point types, trajectory types, and global functions ---

constexpr size_t g_DIMENSIONS = 2;
template <size_t dims> using _point_t = std::array<double, dims>;

struct get_point_coord {
  template <size_t dim> static double get(const _point_t<g_DIMENSIONS> &p) {
    return std::get<dim>(p);
  }
};

template <size_t dims> using _trajectory_t = std::vector<_point_t<dims>>;

template <size_t dims>
void _copyMatrixToTrajectory(const NumericMatrix &m, _trajectory_t<dims> &t) {
  for (size_t i = 0; i < m.nrow(); i++) {
    for (size_t d = 0; d < dims; d++)
      t[i][d] = m(i, d);
  }
}

template <size_t dims>
NumericMatrix _trajectory_to_matrix(const _trajectory_t<dims> &t) {
  NumericMatrix m(t.size(), dims);
  for (size_t i = 0; i < t.size(); i++)
    for (size_t d = 0; d < dims; d++)
      m(i, d) = t[i][d];
  return m;
}

template <size_t dims> class append_to_list {
  List *_list;

public:
  append_to_list(List &l) : _list(&l) {}

  void operator()(const _trajectory_t<dims> &t) {
    _list->push_back(_trajectory_to_matrix<dims>(t));
  }
};

#define ASSERT_VALID_DATASET(dataset, k)                                       \
  {                                                                            \
    if (k < 0 || k >= dataset.size())                                          \
      throw(std::runtime_error("Invalid handle!"));                            \
  }

// ------------------- duetschvahrenhold::grid -------------------

/// grid_data caches the data outside of R to have valid container
/// as some of the implementations are not yet compatible.
template <size_t dims> class grid_data {
public:
  typedef frechetrange::detail::dv::grid<
      dims, _trajectory_t<dims>, get_point_coord>
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
  std::vector<_trajectory_t<dims>> _trajectories;
  grid_type *_pGrid;
};

std::vector<grid_data<g_DIMENSIONS>> g_grids;

// [[Rcpp::export]]
size_t internal_dv_create_index() {
  g_grids.emplace_back();
  return g_grids.size() - 1;
}

// [[Rcpp::export]]
size_t internal_dv_add_trajectory(size_t handle, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_grids, handle);
  return g_grids[handle].add_trajectory(m);
}

// [[Rcpp::export]]
size_t internal_dv_clear(size_t handle) {
  ASSERT_VALID_DATASET(g_grids, handle);
  g_grids[handle].clear();
}

// [[Rcpp::export]]
bool internal_dv_build_index(size_t handle, double meshSize) {
  ASSERT_VALID_DATASET(g_grids, handle);
  Rcout << "Creating dv index from " << g_grids[handle].size()
        << " trajectories." << std::endl;
  g_grids[handle].build_index(meshSize);
}

// [[Rcpp::export]]
List internal_dv_range_query(size_t handle, const NumericMatrix &m,
                             double dist) {
  ASSERT_VALID_DATASET(g_grids, handle);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);

  List resultSet;
  append_to_list<g_DIMENSIONS> append(resultSet);
  g_grids[handle].grid().range_query(t, dist, append);
  return resultSet;
}

// ------------------- baldusbringmann::spatial_index -------------------

template <size_t dims>
using bb_index_type = frechetrange::detail::bb::spatial_index<
    g_DIMENSIONS, _trajectory_t<dims>, get_point_coord>;

std::vector<bb_index_type<g_DIMENSIONS>> g_treeIndices;

// [[Rcpp::export]]
size_t internal_bb_create_index() {
  g_treeIndices.emplace_back();
  return g_treeIndices.size() - 1;
}

// [[Rcpp::export]]
size_t internal_bb_add_trajectory(size_t handle, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_treeIndices, handle);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  g_treeIndices[handle].insert(t);
  return (g_treeIndices[handle].size() - 1);
}

// [[Rcpp::export]]
List internal_bb_range_query(size_t handle, const NumericMatrix &m,
                             double dist) {
  ASSERT_VALID_DATASET(g_treeIndices, handle);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);

  List resultSet;
  append_to_list<g_DIMENSIONS> append(resultSet);
  g_treeIndices[handle].range_query(t, dist, append);
  return resultSet;
}

// ------------------- tue::spatial_hash -------------------

template <size_t dims>
using tue_index_type =
    frechetrange::detail::bddm::spatial_hash<g_DIMENSIONS, _trajectory_t<dims>,
                                            get_point_coord>;

std::vector<tue_index_type<g_DIMENSIONS>> g_tue_indices;

// [[Rcpp::export]]
size_t internal_tue_create_index() {
  g_tue_indices.emplace_back();
  return g_tue_indices.size() - 1;
}

// [[Rcpp::export]]
size_t internal_tue_add_trajectory(size_t handle, const NumericMatrix &m) {
  ASSERT_VALID_DATASET(g_tue_indices, handle);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);
  g_tue_indices[handle].insert(std::move(t));
  return (g_tue_indices[handle].size() - 1);
}

// [[Rcpp::export]]
bool internal_tue_build_index(size_t handle) {
  ASSERT_VALID_DATASET(g_tue_indices, handle);
  Rcout << "Creating tue index from " << g_tue_indices[handle].size()
        << " trajectories." << std::endl;
  g_tue_indices[handle].build_index();
}

// [[Rcpp::export]]
List internal_tue_range_query(size_t handle, const NumericMatrix &m,
                              double dist) {
  ASSERT_VALID_DATASET(g_tue_indices, handle);
  _trajectory_t<g_DIMENSIONS> t(m.nrow());
  _copyMatrixToTrajectory<g_DIMENSIONS>(m, t);

  List resultSet;
  append_to_list<g_DIMENSIONS> append(resultSet);
  g_tue_indices[handle].range_query(t, dist, append);
  return resultSet;
}
