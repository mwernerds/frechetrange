#ifndef BALDUS_BRINGMANN_SPATIAL_INDEX_INC
#define BALDUS_BRINGMANN_SPATIAL_INDEX_INC

// A point in the n-dimensional space
template <size_t dimensions>
using nd_point = std::array<distance_t, dimensions>;

// ----------- quadtree -----------

template <size_t dimensions>
distance_t nd_point_dist(const nd_point<dimensions> &a,
                         const nd_point<dimensions> &b);

template <>
distance_t nd_point_dist(const nd_point<8> &a, const nd_point<8> &b) {
  distance_t result = 0;
  for (size_t i = 0; i < 4; i += 2) {
    result = std::max(result,
                      std::sqrt(sqr(a[i] - b[i]) + sqr(a[i + 1] - b[i + 1])));
  }
  for (size_t i = 4; i < 8; ++i) {
    result = std::max(result, std::abs(a[i] - b[i]));
  }
  return result;
}

constexpr long long constexpr_power(long long base, long long exponent) {
  return (exponent == 0) ? 1 : base * constexpr_power(base, exponent - 1);
}

// Returns true iff the i-th lowest bit (starting at 0) in number is 1.
constexpr bool bit_is_set(long long number, int i) {
  return (number & (1 << i)) != 0;
}

template <size_t dimensions, typename element, size_t max_elments_per_node>
struct quadtree_node {
  typedef nd_point<dimensions> point;

  std::array<std::unique_ptr<quadtree_node>, constexpr_power(2, dimensions)>
      _children;
  point _min_coords, _max_coords;
  std::vector<std::pair<point, element>> _elements;
  bool _has_children;

  quadtree_node(const point &min_coords, const point &max_coords)
      : _children(), _min_coords(min_coords), _max_coords(max_coords),
        _elements(), _has_children(false) {}

  std::vector<element> get(const point &center, distance_t radius) const {
    std::vector<element> result;
    get_recursive(center, radius, result);
    return result;
  }

  void get_recursive(const point &center, distance_t radius,
                     std::vector<element> &result) const {
    bool covered = true;
    for (size_t i = 0; i < dimensions; ++i) {
      if (_min_coords[i] > center[i] + radius ||
          _max_coords[i] < center[i] - radius) {
        return; // circle does not intersect current area
      }
      if (_min_coords[i] < center[i] - radius ||
          _max_coords[i] > center[i] + radius) {
        covered = false;
      }
    }

    if (_has_children && !covered) {
      for (auto &child : _children) {
        child->get_recursive(center, radius, result);
      }

    } else {
      for (const std::pair<point, element> &e : _elements) {
        if (nd_point_dist(center, e.first) <= radius) {
          result.push_back(e.second);
        }
      }
    }
  }

  size_t get_child_id(const point &position) const {
    size_t result = 0;
    for (int i = dimensions - 1; i >= 0; --i) {
      bool high = position[i] >= (_min_coords[i] + _max_coords[i]) / 2;
      result = 2 * result + high;
    }
    return result;
  }

  void add(const element &e, const point &position) {
    _elements.push_back(std::make_pair(position, e));
    if (_has_children) {
      _children[get_child_id(position)]->add(e, position);
    } else {
      if (_elements.size() > max_elments_per_node)
        split();
    }
  }

  void split() {
    bool splittable = false;
    for (size_t i = 1; i < _elements.size(); ++i) {
      if (_elements[i].first != _elements[i - 1].first) {
        splittable = true;
        break;
      }
    }

    if (!splittable)
      return;

    for (size_t i = 0; i < _children.size(); ++i) {
      point mini = _min_coords, maxi = _max_coords;
      for (size_t j = 0; j < dimensions; ++j) {
        distance_t center = (_min_coords[j] + _max_coords[j]) / 2;
        if (bit_is_set(i, j)) {
          mini[j] = center;
        } else {
          maxi[j] = center;
        }
      }
      _children[i] =
          std::unique_ptr<quadtree_node>(new quadtree_node(mini, maxi));
    }

    _has_children = true;

    for (auto &e : _elements) {
      _children[get_child_id(e.first)]->add(e.second, e.first);
    }
  }
};

template <size_t dimensions, typename element,
          size_t max_elments_per_node = constexpr_power(2, dimensions)>
using quadtree = quadtree_node<dimensions, element, max_elments_per_node>;

// ----------- spatial_index -----------

template <typename Trajectory, typename squareddistancefunctional,
          typename xgetterfunctional, typename ygetterfunctional>
class spatial_index {
public:
  spatial_index(squareddistancefunctional dist2, xgetterfunctional xGetter,
                ygetterfunctional yGetter)
      : _q({{min_x, min_y, min_x, min_y, min_x, min_y, min_x, min_y}},
           {{max_x, max_y, max_x, max_y, max_x, max_y, max_x, max_y}}),
        _curves(), _decider(dist2, xGetter, yGetter), _dist2(dist2),
        _getX(xGetter), _getY(yGetter) {}

  void add_curve(const Trajectory &t) {
    _q.add(_curves.size(), get_position(t));
    _curves.emplace_back(t, _dist2);
  }

  // Returns the ids of all trajectories that may have frechet distance d or
  // less.
  // The result may however contain some whose frechet distance to c is too
  // large.
  std::vector<const Trajectory *> get_close_curves(const Trajectory &t,
                                                   distance_t d) const {
    std::vector<const Trajectory *> resultSet;
    auto pushBackResult = [&resultSet](const Trajectory *t) {
      resultSet.push_back(t);
    };
    get_close_curves(t, d, pushBackResult);
    return resultSet;
  }

  template <typename OutputFunctional>
  void get_close_curves(const Trajectory &t, distance_t d,
                        OutputFunctional &output) const {
    // TODO: don't copy t
    curve<Trajectory> c(t, _dist2);
    std::vector<size_t> potential_curves = _q.get(get_position(t), d);

    for (size_t i : potential_curves) {
      const curve<Trajectory> &c2 = _curves[i];

      if (get_frechet_distance_upper_bound(t, c2.trajectory()) <= sqr(d)) {
        output(&(c2.trajectory()));
      } else if (negfilter(c, c2, d)) {
        continue;
      } else if (_decider.is_frechet_distance_at_most(c, c2, d)) {
        output(&(c2.trajectory()));
      }
    }
  }

private:
  static constexpr distance_t min_x = -1e18l, min_y = -1e18l, max_x = 1e18l,
                              max_y = 1e18l;

  quadtree<8, size_t>
      _q; // first x, first y, last x, last y, min x, min y, max x, max y
  std::vector<curve<Trajectory>> _curves;

  FrechetDistance<squareddistancefunctional, xgetterfunctional,
                  ygetterfunctional>
      _decider;
  squareddistancefunctional _dist2;
  xgetterfunctional _getX;
  ygetterfunctional _getY;

  template <typename Point> distance_t dist(Point &p, Point &q) {
    return std::sqrt(_dist2(p, q));
  }

  nd_point<8> get_position(const Trajectory &t) const {
    size_t last = t.size() - 1;
    nd_point<8> p = {{_getX(t[0]), _getY(t[0]), _getX(t[last]), _getY(t[last]),
                      max_x, max_y, min_x, min_y}};
    for (size_t i = 0; i <= last; ++i) {
      p[4] = std::min(p[4], _getX(t[i]));
      p[5] = std::min(p[5], _getY(t[i]));
      p[6] = std::max(p[6], _getX(t[i]));
      p[7] = std::max(p[7], _getY(t[i]));
    }
    return p;
  }

  // ----------- filters -----------

  /*
   * Calculates an upper bound for the squared Frechet distance of a and b by
   * guessing a matching between a and b
   * O(a.size() + b.size())
   */
  distance_t get_frechet_distance_upper_bound(const Trajectory &a,
                                              const Trajectory &b) const {
    distance_t sqrdDist = _dist2(a[a.size() - 1], b[b.size() - 1]);
    size_t pos_a = 0, pos_b = 0;

    while (pos_a + pos_b < a.size() + b.size() - 2) {
      sqrdDist = std::max(sqrdDist, _dist2(a[pos_a], b[pos_b]));
      if (pos_a == a.size() - 1)
        ++pos_b;
      else if (pos_b == b.size() - 1)
        ++pos_a;
      else {
        distance_t dist_a = _dist2(a[pos_a + 1], b[pos_b]);
        distance_t dist_b = _dist2(a[pos_a], b[pos_b + 1]);
        distance_t dist_both = _dist2(a[pos_a + 1], b[pos_b + 1]);
        if (dist_a < dist_b && dist_a < dist_both) {
          ++pos_a;
        } else if (dist_b < dist_both) {
          ++pos_b;
        } else {
          ++pos_a;
          ++pos_b;
        }
      }
    }

    return sqrdDist;
  }

  /*
   * Returns (a discrete approximation of) the first point c[j] on c, with j >=
   * i, that is within distance d of point p.
   */
  template <typename Point>
  size_t nextclosepoint(const curve<Trajectory> &c, size_t i, const Point &p,
                        distance_t d) const {
    const Trajectory &t = c.trajectory();
    size_t delta = 1;
    size_t k = i;
    while (true) {
      if (k == t.size() - 1) {
        if (_dist2(t[k], p) <= sqr(d)) {
          return k;
        } else {
          return c.size();
        }
      } else {
        delta = std::min(delta, t.size() - 1 - k);
        if (std::sqrt(_dist2(p, t[k])) - c.curve_length(k, k + delta) > d) {
          k += delta;
          delta *= 2;
        } else if (delta > 1) {
          delta /= 2;
        } else {
          return k;
        }
      }
    }
  }

  /*
   * Tries to show that the Frechet distance of a and b is more than d. Returns
   * true if a proof is found.
   */
  bool negfilter(const curve<Trajectory> &c1, const curve<Trajectory> &c2,
                 distance_t d) const {
    for (size_t delta = std::max(c1.size(), c2.size()) - 1; delta >= 1;
         delta /= 2) {
      size_t i = 0;
      for (size_t j = 0; j < c2.size(); j += delta) {
        i = nextclosepoint(c1, i, c2.trajectory()[j], d);
        if (i >= c1.size()) {
          return true;
        }
      }
      size_t j = 0;
      for (size_t i = 0; i < c1.size(); i += delta) {
        j = nextclosepoint(c2, j, c1.trajectory()[i], d);
        if (j >= c2.size()) {
          return true;
        }
      }
    }
    return false;
  }
};

#endif
