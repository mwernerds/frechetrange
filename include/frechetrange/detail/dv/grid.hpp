#ifndef DV_GRID_HPP
#define DV_GRID_HPP

#include <algorithm>  // for std::sort, std::lower_bound, std::upper_bound, and std::max
#include <array>
#include <cmath>       // for std::floor
#include <functional>  // for std::hash
#include <stdexcept>   // for std::invalid_argument
#include <unordered_map>
#include <utility>  // for std::move
#include <vector>

#ifdef ENABLE_MULTITHREADING
#include <deque>
#include <future>  // for std::async
#endif

#include "frechet_distance.hpp"

namespace frechetrange {
namespace detail {
namespace dv {

/**
* A grid of fixed mesh size spanning the Euclidean plane.
* It stores dataset trajectories in its cells to
* efficiently perform range queries on their MBR corners.
* Trajectories are mapped to cells using one of the four bounding box corners.
*/
// TODO: multidimensional
template <size_t dimensions, typename trajectory, typename get_coordinate,
          typename squared_distance =
              euclidean_distance_sqr<dimensions, get_coordinate>>
class grid {
   public:
    /**
    * Creates a grid of the specified mesh size.
    */
    grid(double mesh_size, const squared_distance &dist2 = squared_distance())
        : _mesh_size(mesh_size),
          _maps(),
          _expected_query_cost{{0, 0, 0, 0}},
          _use_left_border(),
          _use_bottom_border(),
          _optimized(false),
          _dist2(dist2) {}
    grid(const grid &) = default;
    grid(grid &&) = default;
    grid &operator=(const grid &) = default;
    grid &operator=(grid &&) = default;

    /**
    * Returns the mesh size of this grid.
    */
    double mesh_size() const { return _mesh_size; }

    /**
    * Reserves internal storage so that the indicated number
    * of inserts can be performed without resizing.
    */
    void reserve(size_t num_trajectories) {
        for (size_t i = 0; i < NUM_MAPS; ++i) {
            if (!_optimized ||
                i == to_map_index(_use_left_border, _use_bottom_border))
                _maps[i].reserve(num_trajectories);
        }
    }

    void insert(const trajectory &t) {
        if (t.size() > 0) insert_impl(mbr_t(t));
    }

    void insert(trajectory &&t) {
        if (t.size() > 0) insert_impl(mbr_t(std::move(t)));
    }

    /**
    * Determines which MBR corner to use during hashing to minimize
    * the expected query cost. Furthermore, the grid cells are sorted
    * to achieve better query times.
    */
    void build_index() {
        if (_optimized) return;

        // only keep the grid with the best expected performance
        choose_best_map();
        for (size_t i = 0; i < NUM_MAPS; ++i) {
            if (i != to_map_index(_use_left_border, _use_bottom_border))
                _maps[i].clear();
        }

        // sort the grid cells
        if (_use_left_border)
            if (_use_bottom_border)
                sort_cells<true, true>();
            else
                sort_cells<true, false>();
        else if (_use_bottom_border)
            sort_cells<false, true>();
        else
            sort_cells<false, false>();

        _optimized = true;
    }

    std::vector<const trajectory *> range_query(const trajectory &query,
                                                double dist_threshold) {
        std::vector<const trajectory *> result_set;
        auto push_back_result = [&result_set](const trajectory &t) {
            result_set.push_back(&t);
        };
        range_query(query, dist_threshold, push_back_result);
        return result_set;
    }

    /**
    * @param output Supports the method operator()(const trajectory &)
    *               to output the result trajectories.
    */
    template <typename output_func>
    void range_query(const trajectory &query, double dist_threshold,
                     output_func &output) {
        if (query.size() == 0)
            return;
        else if (dist_threshold > _mesh_size)
            throw std::invalid_argument(
                "The distance threshold is grater than the mesh size.");

        if (!_optimized) choose_best_map();

        if (_use_left_border)
            if (_use_bottom_border)
                this->template query<true, true, output_func>(
                    query, dist_threshold, output);
            else
                this->template query<true, false, output_func>(
                    query, dist_threshold, output);
        else if (_use_bottom_border)
            this->template query<false, true, output_func>(
                query, dist_threshold, output);
        else
            this->template query<false, false, output_func>(
                query, dist_threshold, output);
    }

   private:
    /**
    * Integer coordinates of a grid cell
    */
    struct cell_nr {
        long long x, y;
        bool operator==(const cell_nr &other) const {
            return x == other.x && y == other.y;
        }
    };

    struct mbr_t {
        // TODO: evaluate using pointer instead
        trajectory t;
        /**
        * Coordinates of the borders of the bounding box
        */
        double min_x, max_x, min_y, max_y;

        mbr_t(const trajectory &t) : t(t) { init_borders(); }
        mbr_t(trajectory &&t) : t(std::move(t)) { init_borders(); }
        mbr_t() = default;
        mbr_t(const mbr_t &) = default;
        mbr_t(mbr_t &&) = default;
        mbr_t &operator=(const mbr_t &) = default;
        mbr_t &operator=(mbr_t &&) = default;

        void init_borders() {
            min_x = max_x = get_coordinate::template get<0>(t[0]);
            min_y = max_y = get_coordinate::template get<1>(t[0]);

            for (size_t i = 1; i < t.size(); ++i) {
                auto x = get_coordinate::template get<0>(t[i]);
                if (x < min_x)
                    min_x = x;
                else if (x > max_x)
                    max_x = x;

                auto y = get_coordinate::template get<1>(t[i]);
                if (y < min_y)
                    min_y = y;
                else if (y > max_y)
                    max_y = y;
            }
        }
        /**
        * Returns the coordinate of the specified border of the bounding box.
        * @tparam x_dim Whether a x- or y-coordinate is returned
        * @tparam first Whether the first (i. e., left resp. bottom) or
        *               second (i.e., right resp. top) border is returned
        */
        template <bool x_dim, bool first>
        double get_border() const {
            // the compiler hopefully eliminates the branches
            if (x_dim)
                if (first)
                    return min_x;
                else
                    return max_x;
            else if (first)
                return min_y;
            else
                return max_y;
        }
    };

    /**
    * Compare MBRs by the specified border coordinate.
    */
    template <bool x_dim, bool first>
    struct mbr_comparator {
        bool operator()(const mbr_t &m1, const mbr_t &m2) const {
            return m1.template get_border<x_dim, first>() <
                   m2.template get_border<x_dim, first>();
        }
    };

    /**
    * Cell of the grid
    */
    struct cell_t {
        using bucket_t = std::vector<mbr_t>;
        /**
        * Whether this cell is sorted by x- or y-coordinates
        */
        bool x_sorted;
        /**
        * List of dataset trajectories stored in this grid cell
        */
        bucket_t bucket;

        cell_t() = default;
        cell_t(const cell_t &) = default;
        cell_t(cell_t &&) = default;
        cell_t &operator=(const cell_t &) = default;
        cell_t &operator=(cell_t &&) = default;

        template <bool left, bool bottom>
        void sort() {
            // decide whether to sort by x- or y-coordinates
            x_sorted = choose_sorting_order<left, bottom>();
            if (x_sorted)
                std::sort(bucket.begin(), bucket.end(),
                          mbr_comparator<true, left>());
            else
                std::sort(bucket.begin(), bucket.end(),
                          mbr_comparator<false, bottom>());
        }

        template <bool left, bool bottom>
        bool choose_sorting_order() {
            // find extremal MBR corner coordinates
            double min_x = bucket[0].template get_border<true, left>();
            double max_x = min_x;
            double min_y = bucket[0].template get_border<false, bottom>();
            double max_y = min_y;

            for (size_t i = 1; i < bucket.size(); ++i) {
                double x = bucket[i].template get_border<true, left>();
                if (x < min_x) {
                    min_x = x;
                } else if (x > max_x) {
                    max_x = x;
                }
                double y = bucket[i].template get_border<false, bottom>();
                if (y < min_y) {
                    min_y = y;
                } else if (y > max_y) {
                    max_y = y;
                }
            }

            // choose the dimension with greater diffusion as sorting order
            double x_range = max_x - min_x;
            double y_range = max_y - min_y;
            return x_range >= y_range;
        }
    };

    struct cell_hasher {
        size_t operator()(const cell_nr &cell) const {
            std::hash<long long> long_hasher;
            size_t hash = long_hasher(cell.x) + 0x9e3779b9;
            hash ^=
                long_hasher(cell.y) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            return hash;
        }
    };

    using frechet_decider =
        frechet_distance<dimensions, get_coordinate, squared_distance>;

    /**
    * Minimal number of containing trajectories for a cell to be sorted.
    */
    static constexpr size_t MIN_SORT_SIZE = 16;
    /**
    * Describes which horizontal resp. vertical cell borders are crossed
    * by a range query
    */
    using crossings = unsigned int;
    /**
    * No cell border is crossed
    */
    static constexpr crossings CROSSES_NONE = 0;
    /**
    * The left resp. bottom cell border is crossed
    */
    static constexpr crossings CROSSES_BEGIN = 2;
    /**
    * The right resp. top cell border is crossed
    */
    static constexpr crossings CROSSES_END = 1;

    /**
    * Mesh size of this grid, i. e., the side length of cells
    */
    double _mesh_size;
    using cell_map = std::unordered_map<cell_nr, cell_t, cell_hasher>;
    static constexpr size_t NUM_MAPS = 4;
    /**
    * Map from the 2D coordinates to cells of dataset trajectories
    * for each of the four MBR corners.
    */
    std::array<cell_map, NUM_MAPS> _maps;
    /**
    * Sum of squared cell sizes of each mapping.
    */
    std::array<size_t, NUM_MAPS> _expected_query_cost;
    /**
    * Whether the bounding box corner, according to which queries are mapped to
    * grid cells, lies on the left (or right) border of the bounding box.
    */
    bool _use_left_border;
    /**
    * Whether the bounding box corner, according to which queries are mapped to
    * grid cells, lies on the bottom (or top) border of the bounding box.
    */
    bool _use_bottom_border;
    /**
    * Whether the cells have been sorted.
    */
    bool _optimized;

    squared_distance _dist2;

    long long to_cell_nr(double point_coord) const {
        return static_cast<long long>(std::floor(point_coord / _mesh_size));
    }

    double to_cell_coord(double point_coord) const {
        return std::floor(point_coord / _mesh_size) * _mesh_size;
    }

    size_t to_map_index(bool left, bool bottom) const {
        return 2 * static_cast<size_t>(left) + static_cast<size_t>(bottom);
    }

    void insert_impl(mbr_t &&mbr) {
        if (!_optimized) {
            insert_in_map<false, false>(mbr_t(mbr));
            insert_in_map<false, true>(mbr_t(mbr));
            insert_in_map<true, false>(mbr_t(mbr));
            insert_in_map<true, true>(std::move(mbr));
        } else {
            if (_use_left_border)
                if (_use_bottom_border)
                    insert_in_map<true, true>(std::move(mbr_t(mbr)));
                else
                    insert_in_map<true, false>(std::move(mbr_t(mbr)));
            else if (_use_bottom_border)
                insert_in_map<false, true>(std::move(mbr_t(mbr)));
            else
                insert_in_map<false, false>(std::move(mbr_t(mbr)));
        }
    }

    template <bool left, bool bottom>
    void insert_in_map(mbr_t &&mbr) {
        // integer coordinates of the grid cell
        cell_nr nr{to_cell_nr(mbr.template get_border<true, left>()),
                   to_cell_nr(mbr.template get_border<false, bottom>())};
        cell_t &cell = _maps[to_map_index(left, bottom)][nr];
        typename cell_t::bucket_t &bucket = cell.bucket;

        if (!_optimized) {
            // append the trajectory to the cell's bucket
            bucket.emplace_back(std::move(mbr));
            // update the sum of squared bucket sizes:
            // (b + 1)^2 = b^2 + 2*(b+1) - 1
            _expected_query_cost[to_map_index(left, bottom)] +=
                2 * bucket.size() - 1;
        } else {
            // insert in the sorted bucket
            auto insert_pos =
                cell.x_sorted
                    ? std::upper_bound(bucket.begin(), bucket.end(), mbr,
                                       mbr_comparator<true, left>())
                    : std::upper_bound(bucket.begin(), bucket.end(), mbr,
                                       mbr_comparator<false, bottom>());
            bucket.emplace(insert_pos, std::move(mbr));
        }
    }

    void choose_best_map() {
        // find the map with best expected query cost
        size_t lowest_cost = _expected_query_cost[0];
        _use_left_border = _use_bottom_border = false;

        if (_expected_query_cost[1] < lowest_cost) {
            lowest_cost = _expected_query_cost[1];
            _use_left_border = false;
            _use_bottom_border = true;
        }

        if (_expected_query_cost[2] < lowest_cost) {
            lowest_cost = _expected_query_cost[2];
            _use_left_border = true;
            _use_bottom_border = false;
        }

        if (_expected_query_cost[3] < lowest_cost) {
            lowest_cost = _expected_query_cost[3];
            _use_left_border = true;
            _use_bottom_border = true;
        }
    }

    template <bool left, bool bottom>
    void sort_cells() {
#ifdef ENABLE_MULTITHREADING
        constexpr size_t MAX_NUM_THREADS = 64;
        std::deque<std::future<void>> threads;
#endif

        for (auto &iter : _maps[to_map_index(left, bottom)]) {
            cell_t &cell = iter.second;
            if (cell.bucket.size() >= MIN_SORT_SIZE) {
#ifdef ENABLE_MULTITHREADING
                threads.push_back(
                    std::async(std::launch::async,
                               &cell_t::template sort<left, bottom>, &cell));
                if (threads.size() >= MAX_NUM_THREADS) {
                    threads.pop_front();
                }
#else
                cell.template sort<left, bottom>();
#endif
            }
        }
    }

    template <bool left, bool bottom, typename output>
    void query(const trajectory &query, double threshold, output &out) const {
        mbr_t query_mbr(query);

        // check which horizontal neighbor cells need to be visited
        double cell_coord_x =
            to_cell_coord(query_mbr.template get_border<true, left>());
        bool visit_left =
            query_mbr.template get_border<true, left>() - threshold <
            cell_coord_x;
        bool visit_right =
            query_mbr.template get_border<true, left>() + threshold >=
            cell_coord_x + _mesh_size;

        // check which vertical neighbor cells need to be visited
        double cell_coord_y =
            to_cell_coord(query_mbr.template get_border<false, bottom>());
        bool visit_bottom =
            query_mbr.template get_border<false, bottom>() - threshold <
            cell_coord_y;
        bool visit_top =
            query_mbr.template get_border<false, bottom>() + threshold >=
            cell_coord_y + _mesh_size;

        // memorize the crossed cell borders
        crossings crossed_verticals =
            static_cast<crossings>(visit_left) * CROSSES_BEGIN +
            static_cast<crossings>(visit_right) * CROSSES_END;
        crossings crossed_horizontals =
            static_cast<crossings>(visit_bottom) * CROSSES_BEGIN +
            static_cast<crossings>(visit_top) * CROSSES_END;

        // integral coordinates of the cell the query trajectory is mapped to
        cell_nr cell_no{
            to_cell_nr(query_mbr.template get_border<true, left>()),
            to_cell_nr(query_mbr.template get_border<false, bottom>())};

        frechet_decider fd(_dist2);

        // visit the center cell
        check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                         crossed_verticals, crossed_horizontals,
                                         fd, out);
        // visit the bottom cell
        --cell_no.y;
        if (visit_bottom)
            check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                             crossed_verticals, CROSSES_END, fd,
                                             out);
        // visit the top cell
        cell_no.y += 2;
        if (visit_top)
            check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                             crossed_verticals, CROSSES_BEGIN,
                                             fd, out);

        --cell_no.x;
        if (visit_left) {
            // visit the top-left cell
            if (visit_top)
                check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                                 CROSSES_END, CROSSES_BEGIN, fd,
                                                 out);
            // visit the left cell
            --cell_no.y;
            check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                             CROSSES_END, crossed_horizontals,
                                             fd, out);
            // visit the bottom-left cell
            --cell_no.y;
            if (visit_bottom)
                check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                                 CROSSES_END, CROSSES_END, fd,
                                                 out);
            cell_no.y += 2;
        }

        cell_no.x += 2;
        if (visit_right) {
            // visit the top-right cell
            if (visit_top)
                check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                                 CROSSES_BEGIN, CROSSES_BEGIN,
                                                 fd, out);
            // visit the right cell
            --cell_no.y;
            check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                             CROSSES_BEGIN, crossed_horizontals,
                                             fd, out);
            // visit the bottom-right cell
            --cell_no.y;
            if (visit_bottom)
                check_cell<left, bottom, output>(cell_no, query_mbr, threshold,
                                                 CROSSES_BEGIN, CROSSES_END, fd,
                                                 out);
        }
    }

    template <bool left, bool bottom, typename output>
    void check_cell(const cell_nr &cell_no, const mbr_t &query_mbr,
                    double threshold, crossings crossed_verticals,
                    crossings crossed_horizontals, frechet_decider &fd,
                    output &out) const {
        const cell_map &map = _maps[to_map_index(left, bottom)];
        // ensure that the cell containts some trajectories
        typename cell_map::const_iterator iter = map.find(cell_no);
        if (iter == map.end()) return;

        const cell_t &cell = iter->second;
        // traverse the contained trajectories according to the sorting order
        if (cell.x_sorted) {
            this->template traverse_bucket<true, left, output>(
                cell.bucket, query_mbr, threshold, crossed_verticals, fd, out);
        } else {
            this->template traverse_bucket<false, bottom, output>(
                cell.bucket, query_mbr, threshold, crossed_horizontals, fd,
                out);
        }
    }

    template <bool x_dim, bool first_border, typename output>
    void traverse_bucket(const typename cell_t::bucket_t &bucket,
                         const mbr_t &query_mbr, double threshold,
                         crossings crossed_borders, frechet_decider &fd,
                         output &out) const {
        if (bucket.size() < MIN_SORT_SIZE || !_optimized) {
            // check each trajectory of this cell,
            // as they are not sorted
            for (const mbr_t &t : bucket) {
                check_trajectory<false, false, false>(query_mbr, threshold, t,
                                                      fd, out);
            }
        } else {  // the trajectories are sorted
            // choose the traversing order and beginning depending on
            // which cell borders are crossed
            if (crossed_borders == CROSSES_END) {
                // traverse the trajectories backwards,
                // until the beginning of the active range is reached
                double active_range_begin =
                    query_mbr.template get_border<x_dim, first_border>() -
                    threshold;
                for (size_t i = bucket.size() - 1;
                     bucket[i].template get_border<x_dim, first_border>() >=
                     active_range_begin;
                     --i) {
                    check_trajectory<true, x_dim, first_border>(
                        query_mbr, threshold, bucket[i], fd, out);

                    if (i == 0) {
                        break;
                    }
                }
            } else {
                // traverse the trajectories forwards,
                // until the end of the active range is reached
                auto search_begin = bucket.cbegin();
                auto search_end = bucket.cend();

                // search for the beginning of the active range,
                // if the range query does not span the first,
                // i. e., left or bottom, cell border
                if (crossed_borders == CROSSES_NONE) {
                    double active_range_begin =
                        query_mbr.template get_border<x_dim, first_border>() -
                        threshold;
                    // decide whether to find the beginning of the active range
                    // by binary searching or linear searching
                    if (use_bin_search(bucket.size(), threshold)) {
                        // binary search for the beginning of the active range
                        search_begin = std::lower_bound(
                            search_begin, search_end, active_range_begin,
                            [](const mbr_t &m, double d) {
                                return m.template get_border<x_dim,
                                                             first_border>() <
                                       d;
                            });
                    } else {
                        // linear search for the beginning of the active range
                        while (
                            search_begin != search_end &&
                            search_begin->template get_border<x_dim,
                                                              first_border>() <
                                active_range_begin) {
                            ++search_begin;
                        }
                    }
                }

                // traverse the trajectories forwards,
                // until the end of the active range is reached
                double activeRangeEnd =
                    query_mbr.template get_border<x_dim, first_border>() +
                    threshold;
                while (
                    search_begin != search_end &&
                    search_begin->template get_border<x_dim, first_border>() <=
                        activeRangeEnd) {
                    check_trajectory<true, x_dim, first_border>(
                        query_mbr, threshold, *search_begin, fd, out);
                    ++search_begin;
                }
            }
        }
    }

    template <bool prechecked, bool x_dim, bool first_border,
              typename output_func>
    void check_trajectory(const mbr_t &query_mbr, double threshold,
                          const mbr_t &traj_mbr, frechet_decider &frechet_dist,
                          output_func &out) const {
        // ensure that all bounding box borders are within range
        if (mbrs_within_range<prechecked, x_dim, first_border>(
                query_mbr, threshold, traj_mbr)) {
            // append the dataset trajectory to the output,
            // if it is within Fréchet distance of the query trajectory
            if (sqrd_farthest_mbr_dist(query_mbr, traj_mbr) <=
                    threshold * threshold ||
                frechet_dist.template is_bounded_by<trajectory>(
                    query_mbr.t, traj_mbr.t, threshold)) {
                const trajectory &result = traj_mbr.t;
                out(result);
            }
        }
    }

    template <bool prechecked, bool x_dim, bool first_border>
    bool mbrs_within_range(const mbr_t &query_mbr, double threshold,
                           const mbr_t &traj_mbr) const {
        if (prechecked) {
            // ensure that the three bounding box borders that have not been
            // checked, yet, are within range
            return borders_within_range<!x_dim, !first_border>(
                       query_mbr, threshold, traj_mbr) &&
                   borders_within_range<x_dim, !first_border>(
                       query_mbr, threshold, traj_mbr) &&
                   borders_within_range<!x_dim, !first_border>(
                       query_mbr, threshold, traj_mbr);
        } else {
            // ensure that all bounding box borders are within range
            return borders_within_range<true, true>(query_mbr, threshold,
                                                    traj_mbr) &&
                   borders_within_range<true, false>(query_mbr, threshold,
                                                     traj_mbr) &&
                   borders_within_range<false, true>(query_mbr, threshold,
                                                     traj_mbr) &&
                   borders_within_range<false, false>(query_mbr, threshold,
                                                      traj_mbr);
        }
    }

    template <bool x_dim, bool first_border>
    bool borders_within_range(const mbr_t &query_mbr, double threshold,
                              const mbr_t &traj_mbr) const {
        return query_mbr.template get_border<x_dim, first_border>() -
                       threshold <=
                   traj_mbr.template get_border<x_dim, first_border>() &&
               traj_mbr.template get_border<x_dim, first_border>() <=
                   query_mbr.template get_border<x_dim, first_border>() +
                       threshold;
    }

    double sqrd_farthest_mbr_dist(const mbr_t &query_mbr,
                                  const mbr_t &traj_mbr) const {
        double dx1 = traj_mbr.template get_border<true, false>() -
                     query_mbr.template get_border<true, true>();
        double dx2 = query_mbr.template get_border<true, false>() -
                     traj_mbr.template get_border<true, true>();
        double dy1 = traj_mbr.template get_border<false, false>() -
                     query_mbr.template get_border<false, true>();
        double dy2 = query_mbr.template get_border<false, false>() -
                     traj_mbr.template get_border<false, true>();
        return std::max(dx1 * dx1, dx2 * dx2) + std::max(dy1 * dy1, dy2 * dy2);
    }

    bool use_bin_search(size_t bucket_size, double threshold) const {
        // minimal number of elements to choose binary over linear search
        constexpr size_t MIN_BINARY_SEARCH_SIZE = 32;
        // expected ratio of elements of this bucket to skip;
        // it is not negative, as _mesh_size > 2*threshold
        // otherwise, a cell border would be crossed and
        // this method would not be called
        double expected_ratio_elems_to_skip = 0.5 - threshold / _mesh_size;
        return expected_ratio_elems_to_skip * bucket_size >=
               MIN_BINARY_SEARCH_SIZE;
    }
};

}  // namespace dv
}  // namespace detail
}  // namespace frechetrange

#endif
