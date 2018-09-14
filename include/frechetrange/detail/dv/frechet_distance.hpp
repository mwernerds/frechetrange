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
        : _dist2(dist2), _left_segment_begins() {}

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
        const double bound_sqrd = distance_bound * distance_bound;

        // ensure that the corners of the free space diagram are reachable
        if (_dist2(traj1[0], traj2[0]) > bound_sqrd ||
            _dist2(traj1[traj1.size() - 1], traj2[traj2.size() - 1]) >
                bound_sqrd)
            return false;

        bool first_is_smaller = traj1.size() <= traj2.size();
        const trajectory &smaller_traj = first_is_smaller ? traj1 : traj2;
        const trajectory &bigger_traj = first_is_smaller ? traj2 : traj1;

        if (smaller_traj.size() == 1) {
            return compare_point_to_trajectory(smaller_traj[0], bigger_traj,
                                               bound_sqrd);
#ifdef USE_POSITIVE_FILTER
        } else if (positive_filter(traj1, traj2, bound_sqrd)) {
            return true;
#endif
        } else if (!match_inner_points_monotonously(traj1, traj2, bound_sqrd) ||
                   !match_inner_points_monotonously(traj2, traj1, bound_sqrd)) {
            // There exists no monotone matchting from one trajectory to the
            // other, such that each point is matchted to a point on a
            // segment within distance.
            return false;
        }

        // use the smaller trajectory as first parameter, to consume less memory
        return traverse_free_space_diagram(smaller_traj, bigger_traj,
                                           bound_sqrd);
    }

   private:
    squared_distance _dist2;

    /**
    * Beginnings of reachable parts of the free space segments on
    * the "frontline", i. e., the right segments of the free space cells lastly
    * processed (and accordingly the left segments of the current column).
    * A segment is reachable, if its respective beginning is <= 1.0.
    */
    mutable std::vector<double> _left_segment_begins;

    /**
    * Decides the Fréchet distance problem for a trajectory consisting of only
    * one point.
    */
    template <typename trajectory, typename point_type>
    bool compare_point_to_trajectory(const point_type &p, const trajectory &t,
                                     const double bound_sqrd) const {
        // the first point has already been tested
        for (size_t i = 1; i < t.size(); ++i) {
            if (_dist2(p, t[i]) > bound_sqrd) {
                return false;
            }
        }
        return true;
    }

#ifdef USE_POSITIVE_FILTER
    template <typename trajectory>
    bool positive_filter(const trajectory &traj1, const trajectory &traj2,
                         const double bound_sqrd) const {
        // Positive greedy filter as developed by Baldus and Bringmann in
        // "A fast implementation of near neighbors queries for Fréchet
        // distance",
        // SIGSPATIAL'17
        size_t idx1 = 0, idx2 = 0;
        while (idx1 < traj1.size() - 1 && idx2 < traj2.size() - 1) {
            // Distances of three next pairings
            double dist1 = _dist2(traj1[idx1 + 1], traj2[idx2]);
            double dist2 = _dist2(traj1[idx1], traj2[idx2 + 1]);
            double dist12 = _dist2(traj1[idx1 + 1], traj2[idx2 + 1]);

            // Find the minimal distance
            if (dist12 < dist1 && dist12 < dist2) {
                if (dist12 > bound_sqrd) {
                    return false;
                }

                ++idx1;
                ++idx2;
            } else if (dist1 < dist2) {
                if (dist12 > bound_sqrd) {
                    return false;
                }

                ++idx1;
            } else {
                if (dist12 > bound_sqrd) {
                    return false;
                }

                ++idx2;
            }
        }

        // Advance to the end of the first trajectory, if necessary
        while (idx1 < traj1.size() - 2) {
            ++idx1;
            if (_dist2(traj1[idx1], traj2[idx2]) > bound_sqrd) {
                return false;
            }
        }

        // Advance to the end of the second trajectory, if necessary
        while (idx2 < traj2.size() - 2) {
            ++idx2;
            if (_dist2(traj1[idx1], traj2[idx2]) > bound_sqrd) {
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
    bool match_inner_points_monotonously(const trajectory &points,
                                         const trajectory &segments,
                                         double bound_sqrd) const {
        // the last point has already been tested
        const size_t points_to_match_end = points.size() - 1;
        const size_t num_segments = segments.size() - 1;

        if (points_to_match_end <= 1 || num_segments == 0) {
            // nothing to test
            return true;
        }

        // the first point has already been tested
        size_t point_idx = 1;
        size_t seg_idx = 0;
        double segment_part = 0.0;
        double begin = 0.0, end = 1.0;
        while (true) {
            // search a point on a segment within distance of the current point
            get_line_circle_intersections(
                segments[seg_idx], segments[seg_idx + 1], points[point_idx],
                bound_sqrd, begin, end);

            if (begin <= 1.0 && end >= segment_part) {
                // found a matching point on the current segment
                if (segment_part < begin) {
                    segment_part = begin;
                }

                // advance to the next point
                ++point_idx;
                if (point_idx == points_to_match_end) {
                    // the last point has already been tested
                    return true;
                }
            } else {
                // advance to the next segment and never go back,
                // so that the matching is monotone
                ++seg_idx;
                segment_part = 0.0;

                if (seg_idx == num_segments) {
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
    bool traverse_free_space_diagram(const trajectory &p1, const trajectory &p2,
                                     const double bound_sqrd) const {
        // init the beginnings of reachable parts of the "frontline",
        // i.e., the right segments of the column lastly processed
        // (and accordingly the left segments of the current column)
        clear_frontline(p1.size() - 1);

        // the bottom-left free space corner has been ensured to be reachable
        _left_segment_begins[0] = 0.0;
        // beginning of the reachable part of the bottommost segment of the
        // current
        // column described as a scalar value
        double bottommost_segment_begin = 0.0;
        // bottommost reachable segment of the current column
        size_t bottommost_reachable_row = 0;

        const size_t last_row = p1.size() - 2;
        const size_t last_column = p2.size() - 2;
        // compute the reachable parts of the free space diagram's cells
        // columnwise
        for (size_t col_idx = 0; col_idx < last_column; ++col_idx) {
            // compute whether the current column's bottommost segment is
            // reachable
            if (bottommost_segment_begin == 0.0 &&
                _dist2(p1[0], p2[col_idx]) > bound_sqrd) {
                bottommost_segment_begin = BEGIN_NOT_REACHABLE;
            }

            // beginning of the reachable part of the current bottom segment
            double curr_bottom_seg_begin = bottommost_segment_begin;
            double curr_bottom_seg_end = curr_bottom_seg_begin;

            double curr_right_seg_begin, curr_right_seg_end;
            // traverse this columns' reachable cells
            for (size_t row_idx = bottommost_reachable_row; row_idx < last_row;
                 ++row_idx) {
                // ensure that this cell is reachable
                if (is_segment_reachable(_left_segment_begins[row_idx]) ||
                    is_segment_reachable(curr_bottom_seg_begin)) {
                    // whether the cell's top right corner lies in free space
                    bool is_top_right_free =
                        _dist2(p1[row_idx + 1], p2[col_idx + 1]) <= bound_sqrd;

                    // compute the reachable part of the current cell's right
                    // segment
                    if (is_top_right_free && curr_bottom_seg_end >= 1.0 &&
                        curr_bottom_seg_begin <= 1.0) {
                        // the entire right segment is reachable
                        curr_right_seg_begin = 0.0;
                        curr_right_seg_end = 1.0;
                    } else {
                        get_line_circle_intersections(
                            p1[row_idx], p1[row_idx + 1], p2[col_idx + 1],
                            bound_sqrd, curr_right_seg_begin,
                            curr_right_seg_end);
                        curr_right_seg_begin = get_reachable_begin(
                            curr_right_seg_begin, curr_right_seg_end,
                            _left_segment_begins[row_idx],
                            curr_bottom_seg_begin);
                    }

                    // compute the reachable part of the current cell's top
                    // segment
                    double curr_top_seg_begin, curr_top_seg_end;
                    if (is_top_right_free &&
                        _left_segment_begins[row_idx + 1] <= 0.0) {
                        // the entire top segment is reachable
                        curr_top_seg_begin = 0.0;
                        curr_top_seg_end = 1.0;
                    } else {
                        get_line_circle_intersections(
                            p2[col_idx], p2[col_idx + 1], p1[row_idx + 1],
                            bound_sqrd, curr_top_seg_begin, curr_top_seg_end);
                        curr_top_seg_begin = get_reachable_begin(
                            curr_top_seg_begin, curr_top_seg_end,
                            curr_bottom_seg_begin,
                            _left_segment_begins[row_idx]);
                    }

                    // add the current cell to the "frontline"
                    _left_segment_begins[row_idx] = curr_right_seg_begin;
                    curr_bottom_seg_begin = curr_top_seg_begin;
                    curr_bottom_seg_end = curr_top_seg_end;
                }

                // update the bottommost reachable segment, if necessary
                if (bottommost_reachable_row == row_idx &&
                    !is_segment_reachable(_left_segment_begins[row_idx])) {
                    ++bottommost_reachable_row;
                }
            }

            // compute the reachable part of the topmost cell's right segment
            if (is_segment_reachable(_left_segment_begins[last_row]) ||
                is_segment_reachable(curr_bottom_seg_begin)) {
                if (_dist2(p1[last_row + 1], p2[col_idx + 1]) <= bound_sqrd &&
                    curr_bottom_seg_end >= 1.0 &&
                    curr_bottom_seg_begin <= 1.0) {
                    // the entire segment is reachable
                    curr_right_seg_begin = 0.0;
                } else {
                    get_line_circle_intersections(
                        p1[last_row], p1[last_row + 1], p2[col_idx + 1],
                        bound_sqrd, curr_right_seg_begin, curr_right_seg_end);
                    curr_right_seg_begin = get_reachable_begin(
                        curr_right_seg_begin, curr_right_seg_end,
                        _left_segment_begins[last_row], curr_bottom_seg_begin);
                }
                _left_segment_begins[last_row] = curr_right_seg_begin;
            }

            // ensure that a segment of the frontline is reachable
            if (bottommost_reachable_row == last_row &&
                !is_segment_reachable(_left_segment_begins[last_row])) {
                return false;
            }
        }

        if (is_segment_reachable(_left_segment_begins[last_row])) {
            // the top-right corner is reachable via its left segment
            return true;
        }

        // traverse the rightmost column
        double curr_bottom_seg_begin = BEGIN_NOT_REACHABLE;
        for (size_t row_idx = bottommost_reachable_row; row_idx < last_row;
             ++row_idx) {
            // ensure that this cell is reachable
            if (is_segment_reachable(_left_segment_begins[row_idx]) ||
                is_segment_reachable(curr_bottom_seg_begin)) {
                // compute the reachable part of the current cell's top segment
                double curr_top_seg_begin, curr_top_seg_end;
                if (_dist2(p1[row_idx + 1], p2[last_column + 1]) <=
                        bound_sqrd &&
                    _left_segment_begins[row_idx + 1] <= 0.0) {
                    // the entire top segment is reachable
                    curr_bottom_seg_begin = 0.0;
                } else {
                    get_line_circle_intersections(
                        p2[last_column], p2[last_column + 1], p1[row_idx + 1],
                        bound_sqrd, curr_top_seg_begin, curr_top_seg_end);
                    curr_bottom_seg_begin = get_reachable_begin(
                        curr_top_seg_begin, curr_top_seg_end,
                        curr_bottom_seg_begin, _left_segment_begins[row_idx]);
                }
            }
        }

        // return whether the top-right corner is reachable via its bottom
        // segment
        return is_segment_reachable(curr_bottom_seg_begin);
    }

    /**
    * Increases the allocated memory of the frontline, if necessary,
    * and marks its segments as not reachable.
    */
    void clear_frontline(size_t rows) const {
        // reserve memory, if necessary
        if (_left_segment_begins.size() < rows) {
            _left_segment_begins.resize(rows, BEGIN_NOT_REACHABLE);
        }

        // mark the frontline as not reachable
        std::fill_n(_left_segment_begins.begin(), rows, BEGIN_NOT_REACHABLE);
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
    void get_line_circle_intersections(const point &p1, const point &p2,
                                       const point &cp,
                                       const double radius_squared,
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
        double c = -radius_squared;
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
            if (_dist2(p1, cp) <= radius_squared) {
                begin = 0.0;
                end = 1.0;
            } else {
                begin = end = BEGIN_NOT_REACHABLE;
            }
            return;
        }

        double sqrt_d = std::sqrt(discriminant);
        begin = (-b - sqrt_d) / (2.0 * a);
        end = (-b + sqrt_d) / (2.0 * a);
    }

    /**
    * Computes beginning of the reachable part of a free space segment.
    * @param curr_segment_begin The beginning of the intersection of free space
    *                         with the segment described as a scalar value.
    * @param curr_segment_end The ending of the intersection of free space
    *                       with the segment described as a scalar value.
    * @param prev_parallel_seg_begin The beginning of the reachable part of the
    *                             parallel segment of the same free space cell
    *                             described as a scalar value.
    * @param prev_orthogonal_seg_begin The beginning of the reachable part of
    * the
    *                               orthogonal segment of the same free space
    *                               cell described as a scalar value.
    * @return A scalar value describing the beginning of the reachable part of
    *         the segment.
    */
    double get_reachable_begin(double curr_segment_begin,
                               double curr_segment_end,
                               double prev_parallel_seg_begin,
                               double prev_orthogonal_seg_begin) const {
        if (curr_segment_end < 0.0) {
            // The segment does not intersect the free space.
            return BEGIN_NOT_REACHABLE;
        } else if (curr_segment_begin > 1.0) {
            // The segment does not intersect free space,
            // and there is no need to change curr_segment_begin.
        } else if (is_segment_reachable(prev_orthogonal_seg_begin)) {
            // The reachable part equals the intersection with the free space.
        } else if (curr_segment_begin < prev_parallel_seg_begin) {
            // The beginning of the reachable part is greater
            // than the beginning of the intersection.
            return (curr_segment_end >= prev_parallel_seg_begin)
                       ? prev_parallel_seg_begin
                       : BEGIN_NOT_REACHABLE;
        }

        return curr_segment_begin;
    }

    /**
    * Returns whether a free space segment with the passed
    * beginning of the reachable part (as returned by get_reachable_begin)
    * is reachable.
    * @param begin The beginning of the reachable part of the segment
    *              described as a scalar value.
    */
    bool is_segment_reachable(double begin) const { return begin <= 1.0; }
};

}  // namespace dv
}  // namespace detail
}  // namespace frechetrange

#endif
