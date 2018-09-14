#ifndef BB_FRECHET_DISTANCE_HPP
#define BB_FRECHET_DISTANCE_HPP

#include <algorithm>  // for std::min, std::max, and std::upper_bound
#include <cassert>
#include <cmath>  // for std::sqrt
#include <limits>
#include <tuple>    // for std::tie and std::ignore
#include <utility>  // for std::pair
#include <vector>

#include "../../distance_sqr.hpp"

namespace frechetrange {
namespace detail {
namespace bb {

typedef double distance_t;

template <typename valtype>
valtype sqr(valtype a) {
    return a * a;
}

/*
 * Represents a trajectory. Additionally to the points given in the input file,
 * we also store the length of any prefix of the trajectory.
 */
template <typename trajectory>
class curve {
    const trajectory _trajectory;
    std::vector<distance_t> _prefix_length;

   public:
    template <typename squared_distance>
    curve(const trajectory &t, const squared_distance &dist2)
        : _trajectory(t), _prefix_length(t.size()) {
        _prefix_length[0] = 0;
        for (size_t i = 1; i < t.size(); ++i)
            _prefix_length[i] =
                _prefix_length[i - 1] + std::sqrt(dist2(t[i - 1], t[i]));
    }

    curve() = default;

    size_t size() const { return _trajectory.size(); }

    const trajectory &get_trajectory() const { return _trajectory; }

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
    typedef std::pair<distance_t, distance_t> interval;  // .first is the
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
    interval intersection_interval(const point &circle_center,
                                   distance_t radius, const point &line_start,
                                   const point &line_end) const {
        // Find points p = line_start + lambda * v with
        //     dist(p, circle_center) = radius
        // <=> sqrt(p.x^2 + p.y^2) = radius
        // <=> p.x^2 + p.y^2 = radius^2
        // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2
        // =
        // radius^2
        // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 *
        // v.x^2)
        // + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 *
        // v.y^2) =
        // radius^2
        // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2
        // *
        // line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
        // let a := v.x^2 + v.y^2, b := 2 * line_start.x * v.x + 2 *
        // line_start.y *
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
            return empty_interval;  // no intersection;
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
    template <typename trajectory, typename point>
    distance_t get_dist_to_point_sqr(const trajectory &t,
                                     const point &p) const {
        distance_t result = 0;
        for (size_t i = 0; i < t.size(); ++i)
            result = std::max(result, _dist2(t[i], p));
        return result;
    }

    template <typename trajectory>
    interval get_reachable_a(size_t i, size_t j, const trajectory &a,
                             const trajectory &b, distance_t d) const {
        distance_t start, end;
        std::tie(start, end) = intersection_interval(a[i], d, b[j], b[j + 1]);
        return {start + j, end + j};
    }

    template <typename trajectory>
    interval get_reachable_b(size_t i, size_t j, const trajectory &a,
                             const trajectory &b, distance_t d) const {
        return get_reachable_a(j, i, b, a, d);
    }

    void merge(std::vector<interval> &v, interval i) const {
        if (is_empty_interval(i)) return;
        if (v.size() && i.first - eps <= v.back().second)
            v.back().second = i.second;
        else
            v.push_back(i);
    }

    template <typename trajectory>
    distance_t get_last_reachable_point_from_start(const trajectory &a,
                                                   const trajectory &b,
                                                   const distance_t d) const {
        size_t j = 0;
        while (j < b.size() - 2 && _dist2(a[0], b[j + 1]) <= sqr(d)) ++j;

        distance_t result;
        tie(std::ignore, result) = get_reachable_a(0, j, a, b, d);
        return result;
    }

    template <typename trajectory>
    void get_reachable_intervals(size_t i_min, size_t i_max, size_t j_min,
                                 size_t j_max, const curve<trajectory> &a,
                                 const curve<trajectory> &b, distance_t d,
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

        if (is_empty_interval(tb) && is_empty_interval(ta)) return;

        const trajectory &t1 = a.get_trajectory();
        const trajectory &t2 = b.get_trajectory();
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
                get_reachable_intervals(i_min, i_max, j_min, split_position, a,
                                        b, d, rb, ra, rb_out, ra_middle);
                get_reachable_intervals(i_min, i_max, split_position, j_max, a,
                                        b, d, rb, ra_middle, rb_out, ra_out);
            } else {
                std::vector<interval> rb_middle;
                size_t split_position = (i_max + i_min) / 2;
                get_reachable_intervals(i_min, split_position, j_min, j_max, a,
                                        b, d, rb, ra, rb_middle, ra_out);
                get_reachable_intervals(split_position, i_max, j_min, j_max, a,
                                        b, d, rb_middle, ra, rb_out, ra_out);
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

    template <typename trajectory>
    bool is_bounded_by(const trajectory &a, const trajectory &b,
                       distance_t d) const {
        return is_bounded_by(curve<trajectory>(a, _dist2),
                             curve<trajectory>(b, _dist2), d);
    }

    template <typename trajectory>
    bool is_bounded_by(const curve<trajectory> &c1, const curve<trajectory> &c2,
                       distance_t d) const {
        assert(c1.size());
        assert(c2.size());

        const trajectory &t1 = c1.get_trajectory();
        const trajectory &t2 = c2.get_trajectory();
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

        get_reachable_intervals(0, c1.size() - 1, 0, c2.size() - 1, c1, c2, d,
                                ra, rb, ra_out, rb_out);

        return ra_out.size() && (ra_out.back().second >=
                                 c2.size() - static_cast<distance_t>(1.5));
    }
};

}  // namespace bb
}  // namespace detail
}  // namespace frechetrange

#endif
