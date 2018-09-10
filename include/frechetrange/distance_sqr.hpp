#ifndef DISTANCE_SQR_HPP
#define DISTANCE_SQR_HPP

#include <cstddef>  // for std::size_t

namespace frechetrange {

namespace detail {

// template meta programming
template <typename point, typename get_coordinate, size_t dim>
struct sum_of_sqr_diffs {
    static double dist(const point &p, const point &q) {
        auto delta_dim = get_coordinate::template get<dim - 1>(p) -
                         get_coordinate::template get<dim - 1>(q);
        return static_cast<double>(delta_dim * delta_dim) +
               sum_of_sqr_diffs<point, get_coordinate, dim - 1>::dist(p, q);
    }
};

template <typename point, typename get_coordinate>
struct sum_of_sqr_diffs<point, get_coordinate, 0> {
    static double dist(const point &, const point &) { return 0.0; }
};

template <typename point, typename get_coordinate, size_t dim>
struct add_sum_of_pairwise_diff_prods {
    static void calc(const point &p, const point &q, const point &m, double &a,
                     double &b, double &c) {
        // u := q - p
        auto u_dim = get_coordinate::template get<dim - 1>(q) -
                     get_coordinate::template get<dim - 1>(p);
        // v := p - m
        auto v_dim = get_coordinate::template get<dim - 1>(p) -
                     get_coordinate::template get<dim - 1>(m);

        a += u_dim * u_dim;
        b += u_dim * v_dim;
        c += v_dim * v_dim;

        add_sum_of_pairwise_diff_prods<point, get_coordinate, dim - 1>::calc(
            p, q, m, a, b, c);
    }
};

template <typename point, typename get_coordinate>
struct add_sum_of_pairwise_diff_prods<point, get_coordinate, 0> {
    static void calc(const point &, const point &, const point &, double &a,
                     double &b, double &c) {}
};

}  // namespace frechetrange

template <size_t dimensions, typename get_coordinate>
struct euclidean_distance_sqr {
    template <typename point>
    double operator()(const point &p, const point &q) const {
        return detail::sum_of_sqr_diffs<point, get_coordinate,
                                        dimensions>::dist(p, q);
    }
};

}  // namespace frechetrange

#endif
