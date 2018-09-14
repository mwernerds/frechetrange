#ifndef BDDM_SPATIAL_HASH_HPP
#define BDDM_SPATIAL_HASH_HPP

#include <algorithm>   // for std::max and std::min
#include <cmath>       // for std::sqrt and std::isfinite
#include <cstddef>     // for std::size_t
#include <functional>  // for std::function
#include <limits>
#include <map>
#include <utility>  // for std::move
#include <vector>

#ifdef ENABLE_MULTITHREADING
#include <thread>
#endif

namespace frechetrange {
namespace detail {
namespace bddm {

/*
    This is the collected and reordered (but not redactionally completed)
implementation of

1. Tom van Diggelen    t.w.t.v.diggelen@student.tue.nl
2. Yago Diez           contact@yagodiez.com
3. Kevin Buchin        k.a.buchin@tue.nl
4. Wouter Meulemans    w.meulemans@tue.nl

submitted to the ACM SIGSPATIAL GIS Cup 2017

  Refactored by Fabian and Martin
*/
// TODO: distance functional
// TODO: multidimensional
template <typename point, size_t dimensions, typename trajectory,
          typename get_coordinate,
          typename squared_distance =
              euclidean_distance_sqr<dimensions, get_coordinate>>
class spatial_hash {
    typedef std::vector<point> vertices;
    struct bounding_box_t;
    struct portal;
    struct ext_trajectory;
    struct simplification;

    class di_hash;
    class agarwal_simplification;
    class progressive_agarwal;
    class cfdq_shortcuts;

   public:
    spatial_hash(const squared_distance &dist2 = squared_distance())
        : _dist2(dist2),
          _src_trajectories(),
          _ext_trajectories(),
          _bounding_box(),
          _di_hash(SLOTS_PER_DIMENSION, TOLERANCE),
          _count(0),
          _avgs_bb_ratio{0.0, 0.0, 0.0, 0.0} {}

    size_t size() const { return _src_trajectories.size(); }

    void insert(const trajectory &t) {
        if (t.size() > 0) {
            _src_trajectories.push_back(t);
            _ext_trajectories.emplace_back(t, _ext_trajectories.size(), _dist2);
        }
    }

    void insert(trajectory &&t) {
        if (t.size() > 0) {
            _ext_trajectories.emplace_back(t, _ext_trajectories.size(), _dist2);
            _src_trajectories.emplace_back(std::move(t));
        }
    }

    // Does all needed preprocessing for the given dataset
    void build_index() {
        construct_simplifications();
        add_pts_to_di_hash();
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
        ext_trajectory query_trajectory(query, 0, _dist2);

        progressive_agarwal agarwal_prog;
        double diagonal = query_trajectory.bounding_box.get_diagonal();
        make_source_simplifications(query_trajectory, query_trajectory,
                                    diagonal, agarwal_prog,
                                    NUM_SIMPLIFICATIONS);
        for (int i = 1; i < NUM_SIMPLIFICATIONS; i++) {
            make_source_simplifications(query_trajectory.simplifications[i],
                                        query_trajectory, diagonal,
                                        agarwal_prog, i - 1);
        }

        cfdq_shortcuts cdfqs;
        collect_di_hash_points(
            query_trajectory, dist_threshold, _dist2,
            [&](const ext_trajectory &t) {
                prune_with_simplifications(
                    cdfqs, query_trajectory, dist_threshold, t,
                    [&](const ext_trajectory &t) {
                        prune_with_equal_time(
                            query_trajectory, dist_threshold, t,
                            [&](const ext_trajectory &t) {
                                prune_with_decision_frechet(
                                    cdfqs, query_trajectory, dist_threshold, t,
                                    _dist2, output);
                            },
                            output);
                    },
                    output);
            });
    }

   private:
    static constexpr int NUM_SIMPLIFICATIONS = 4;
    static constexpr int SLOTS_PER_DIMENSION = 500;
    static constexpr double TOLERANCE = 0.00001;

    squared_distance _dist2;

    std::vector<trajectory> _src_trajectories;
    std::vector<ext_trajectory> _ext_trajectories;
    bounding_box_t _bounding_box;
    di_hash _di_hash;

    int _count;
    double _avgs_bb_ratio[4];

    /* DATA STRUCTURES */

    // Represents a boundingbox around a trajectory or group of trajectories
    struct bounding_box_t {
        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::max();

        double maxx = std::numeric_limits<double>::min();
        double maxy = std::numeric_limits<double>::min();

        void add_point(const point &p) {
            add_point(get_coordinate::template get<0>(p),
                      get_coordinate::template get<1>(p));
        }

        void add_point(double x, double y) {
            if (x < minx) minx = x;
            if (y < miny) miny = y;

            if (x > maxx) maxx = x;
            if (y > maxy) maxy = y;
        }

        double get_diagonal() {
            double width = maxx - minx;
            double height = maxy - miny;
            return std::sqrt(width * width + height * height);
        }
    };

    struct portal {
        int source;
        int destination;
        double distance;

        bool operator<(const portal &p) const {
            return destination < p.destination;
        };
    };

    struct ext_trajectory {
        size_t index;
        vertices points;
        std::vector<double> distances;  // distance between points
        std::vector<double> totals;     // total length at vertex x
        std::vector<size_t>
            sourceIndex;  // mapping back to the source trajectory if
                          // this is a simplification
        bounding_box_t bounding_box;

        std::map<int, std::vector<portal>> simp_portals;  // freespace jumps
        std::vector<simplification> simplifications;

        ext_trajectory() = default;
        ext_trajectory(ext_trajectory &&) = default;

        ext_trajectory(const trajectory &t, size_t idx, squared_distance &dist2)
            : index(idx) {
            points.reserve(t.size());
            points.push_back(t[0]);
            distances.push_back(0.0);
            totals.push_back(0.0);
            sourceIndex.push_back(0);
            bounding_box.add_point(t[0]);

            for (size_t i = 1; i < t.size(); ++i) {
                const point &p = t[i];
                double x = get_coordinate::template get<0>(p);
                double y = get_coordinate::template get<1>(p);

                // ignore NaN and duplicates
                const point &prev = points.back();
                if (std::isfinite(x) && std::isfinite(y) &&
                    (x != get_coordinate::template get<0>(prev) ||
                     y != get_coordinate::template get<1>(prev))) {
                    points.push_back(p);
                    distances.push_back(std::sqrt(dist2(prev, points.back())));
                    totals.push_back(totals.back() + distances.back());
                    sourceIndex.push_back(sourceIndex.size());
                    bounding_box.add_point(p);
                }
            }
        }

        size_t size() const { return points.size(); }

        const point &front() const { return points[0]; }

        const point &back() const { return points.back(); }

        void add_simplification(simplification &&simp) {
            simplifications.emplace_back(std::move(simp));
        }
    };

    struct simplification : public ext_trajectory {
        double simplification_eps;
        std::vector<portal> portals;

        simplification() = default;
        simplification(simplification &&) = default;
        simplification(double eps) : simplification_eps(eps), portals() {}
    };

    static double equal_time_distance(const ext_trajectory &p,
                                      const ext_trajectory &q) {
        return equal_time_distance(p.points, q.points, p.totals, q.totals,
                                   p.distances, q.distances, p.size(), q.size(),
                                   0, 0);
    }

    // Implements Equal Time Distance algorithm between two trajectories.
    // The ETD algorithm computes an approximation of frechet distance by
    // taking the 'dog leash' length when traversing two trajectories at
    // the same speed. Used by agarwal and simplification step.
    // If found relevant, this function can be optimized using SIMD instructions
    static double equal_time_distance(const vertices &pverts,
                                      const vertices &qverts,
                                      const std::vector<double> &ptotals,
                                      const std::vector<double> &qtotals,
                                      const std::vector<double> &pdistances,
                                      const std::vector<double> &qdistances,
                                      size_t psize, size_t qsize, size_t pstart,
                                      size_t qstart) {
        double pdist_offset = ptotals[pstart];
        double qdist_offset = qtotals[qstart];
        double pdist = ptotals[psize - 1] - pdist_offset;
        double qdist = qtotals[qsize - 1] - qdist_offset;
        double pscale = qdist / pdist;

        // startpoints
        double p_ptx = get_coordinate::template get<0>(pverts[pstart]);
        double p_pty = get_coordinate::template get<1>(pverts[pstart]);
        double q_ptx = get_coordinate::template get<0>(qverts[qstart]);
        double q_pty = get_coordinate::template get<1>(qverts[qstart]);
        double smax = (p_ptx - q_ptx) * (p_ptx - q_ptx) +
                      (p_pty - q_pty) * (p_pty - q_pty);
        p_ptx = get_coordinate::template get<0>(pverts[psize - 1]);
        p_pty = get_coordinate::template get<1>(pverts[psize - 1]);
        q_ptx = get_coordinate::template get<0>(qverts[qsize - 1]);
        q_pty = get_coordinate::template get<1>(qverts[qsize - 1]);
        double emax = (p_ptx - q_ptx) * (p_ptx - q_ptx) +
                      (p_pty - q_pty) * (p_pty - q_pty);
        if (qdist == 0 || pdist == 0) return std::sqrt(std::max(emax, smax));

        double position = 0;  // from 0 to 1
        size_t p_ptr = pstart + 1;
        size_t q_ptr = qstart + 1;

        while (!(p_ptr == psize - 1 && q_ptr == qsize - 1)) {
            double posP = position * pdist;
            double posQ = position * qdist;
            double nextDistP = ptotals[p_ptr] - pdist_offset - posP;
            double nextDistQ = qtotals[q_ptr] - qdist_offset - posQ;

            if (p_ptr == psize - 1)
                nextDistP = std::numeric_limits<double>::max();
            if (q_ptr == qsize - 1)
                nextDistQ = std::numeric_limits<double>::max();

            if (nextDistP * pscale < nextDistQ) {  // treat P first
                p_ptx = get_coordinate::template get<0>(pverts[p_ptr]);
                p_pty = get_coordinate::template get<1>(pverts[p_ptr]);
                position = (ptotals[p_ptr] - pdist_offset) / pdist;
                double scale =
                    (position * qdist - (qtotals[q_ptr - 1] - qdist_offset)) /
                    qdistances[q_ptr];
                double dx = get_coordinate::template get<0>(qverts[q_ptr]) -
                            get_coordinate::template get<0>(qverts[q_ptr - 1]);
                double dy = get_coordinate::template get<1>(qverts[q_ptr]) -
                            get_coordinate::template get<1>(qverts[q_ptr - 1]);
                q_ptx = get_coordinate::template get<0>(qverts[q_ptr - 1]) +
                        dx * scale;
                q_pty = get_coordinate::template get<1>(qverts[q_ptr - 1]) +
                        dy * scale;
                p_ptr++;
            } else {  // treat Q first
                q_ptx = get_coordinate::template get<0>(qverts[q_ptr]);
                q_pty = get_coordinate::template get<1>(qverts[q_ptr]);

                position = (qtotals[q_ptr] - qdist_offset) / qdist;
                double scale =
                    (position * pdist - (ptotals[p_ptr - 1] - pdist_offset)) /
                    pdistances[p_ptr];
                double dx = get_coordinate::template get<0>(pverts[p_ptr]) -
                            get_coordinate::template get<0>(pverts[p_ptr - 1]);
                double dy = get_coordinate::template get<1>(pverts[p_ptr]) -
                            get_coordinate::template get<1>(pverts[p_ptr - 1]);
                p_ptx = get_coordinate::template get<0>(pverts[p_ptr - 1]) +
                        dx * scale;
                p_pty = get_coordinate::template get<1>(pverts[p_ptr - 1]) +
                        dy * scale;
                q_ptr++;
            }
            double nm = (p_ptx - q_ptx) * (p_ptx - q_ptx) +
                        (p_pty - q_pty) * (p_pty - q_pty);
            if (nm > smax) {
                smax = nm;
            }
        }

        // endpoints
        p_ptx = get_coordinate::template get<0>(pverts[p_ptr]);
        p_pty = get_coordinate::template get<1>(pverts[p_ptr]);
        q_ptx = get_coordinate::template get<0>(qverts[q_ptr]);
        q_pty = get_coordinate::template get<1>(qverts[q_ptr]);
        double nm = (p_ptx - q_ptx) * (p_ptx - q_ptx) +
                    (p_pty - q_pty) * (p_pty - q_pty);
        if (nm > smax) {
            smax = nm;
        }

        return std::sqrt(smax);
    }

    struct range {
        double start = 0.0;
        double end = 0.0;

        bool is_complete() const { return start == 0.0 && end == 1.0; }
    };

    static bool compute_interval(const point &a, const point &b1,
                                 const point &b2, double eps, range &r) {
        // compute the interval along b1-b2 that is within eps-range of a:
        // L(t) = b1 + t * (b2 - b1)
        // D(t) = |L(t) - a|
        // for which t, does D(t) = eps hold?
        // square it to make it easier
        // (b1_x + t * (b2_x - b1_x) - a_x)^2
        //    + (b1_y + t * (b2_y - b1_y) - a_y)^2 = eps^2
        // so, we get:
        // A * t^2 + B * t + C = 0 with
        // A = (b2_x - b1_x)^2 + (b2_y - b1_y)^2
        // B = 2 * ((b2_x - b1_x)*(b1_x - a_x) + (b2_y - b1_y)*(b1_y - a_y));
        // C = (b1_x - a_x)^2 + (b1_y - a_y)^2 - eps^2
        double b2m1x = get_coordinate::template get<0>(b2) -
                       get_coordinate::template get<0>(b1);
        double b2m1y = get_coordinate::template get<1>(b2) -
                       get_coordinate::template get<1>(b1);
        double b1max = get_coordinate::template get<0>(b1) -
                       get_coordinate::template get<0>(a);
        double b1may = get_coordinate::template get<1>(b1) -
                       get_coordinate::template get<1>(a);

        double A = b2m1x * b2m1x + b2m1y * b2m1y;
        double B = 2 * ((b2m1x) * (b1max) + (b2m1y) * (b1may));
        double C = b1max * b1max + b1may * b1may - eps * eps;

        double D = B * B - 4 * A * C;
        if (D < 0) {
            // no solution
            return false;
        } else {
            // pull out the Sqrt(D)
            double sqrtD = std::sqrt(D);
            double t1 = (-B + sqrtD) / (2 * A);
            double t2 = (-B - sqrtD) / (2 * A);
            // rather than doing a swap, if may be faster to check the sign of A
            // before doing the assignment, OR just use min/max, etc
            double tempt1 = t1;
            t1 = std::min(t1, t2);
            t2 = std::max(tempt1, t2);

            if (t2 < 0 || t1 > 1) {
                return false;
            } else {
                r.start = std::max(0.0, t1);
                r.end = std::min(1.0, t2);
                return true;
            }
        }
    };

    // Adapted from implementation of Yago Diez, could be sped up but is not a
    // bottleneck
    class di_hash {
        struct vertex {
            const point &p;
            const ext_trajectory &t;
            bool is_start;
        };

        const int _slots_per_dim;  // number of equally spaced subdivisions in
                                   // each dimension
        std::vector<std::vector<vertex>> _elements;
        double
            _tolerance;  // tolerance to prevent numerical representation errors
        double
            _limits[2]
                   [2];  // two rows (x,y) and two columns(min,max), keeps the
                         // information on limits in each dimension, for 3D,
                         // just add one row.

       public:
        di_hash(int numC, double tol)
            : _slots_per_dim(numC), _elements(numC * numC), _tolerance(tol) {}

        void init(const bounding_box_t &bounding_box) {
            // First, locate numerichal limits
            // initialize limits at a three-component two doubles (max and min)
            // vector
            // set initial values for extreme values
            _limits[0][0] = bounding_box.minx;  // min in X
            _limits[0][1] = bounding_box.maxx;  // max in X
            _limits[1][0] = bounding_box.miny;  // min in Y
            _limits[1][1] = bounding_box.maxy;  // max in Y
        }

        void add_point(const point &p, const ext_trajectory &t, bool is_start) {
            int xSlot =
                find_slot(get_coordinate::template get<0>(p), 'x', false);
            int ySlot =
                find_slot(get_coordinate::template get<1>(p), 'y', false);
            _elements[xSlot * _slots_per_dim + ySlot].push_back(
                vertex{p, t, is_start});
        }

        // Return neighbors at range strictly less than distance for a given
        // point,
        // does not return the query point if present
        void neighbors_with_callback(
            const point &start, const point &end, double eps,
            std::vector<ext_trajectory> &trajectories, squared_distance &dist2,
            const std::function<void(const ext_trajectory &)> &emit) {
            auto limitsX =
                slots_touched(get_coordinate::template get<0>(start), eps, 'x');
            auto limitsY =
                slots_touched(get_coordinate::template get<1>(start), eps, 'y');
            for (int i = limitsX.first; i <= limitsX.second; i++) {
                for (int j = limitsY.first; j <= limitsY.second; j++) {
                    for (const vertex &pActual :
                         _elements[i * _slots_per_dim + j]) {
                        if (pActual.is_start) {
                            double distSQ = dist2(start, pActual.p);
                            if (distSQ < eps * eps) {
                                const ext_trajectory &t = pActual.t;
                                if (dist2(end, t.back()) < eps * eps) {
                                    emit(t);
                                }
                            }
                        }
                    }
                }
            }
        }

       private:
        // Given a type of search (x or y) and a search range (min, max(, return
        // the
        // two indexes of the first and last slots affected by the search
        std::pair<int, int> slots_touched(double center, double eps,
                                          char type) {
            double min = center - eps;
            double max = center + eps;
            // first position of the following vector is for the minimum slot
            // affected
            // and second for the maximum
            return std::make_pair(find_slot(min, type, true),
                                  find_slot(max, type, true));
        }

        int find_slot(double val, char type, bool allow_overflow) {
            double min, max;

            switch (type) {
                case 'x':
                    min = _limits[0][0];
                    max = _limits[0][1];
                    break;
                default:  // case 'y':
                    min = _limits[1][0];
                    max = _limits[1][1];
                    break;
            }

            if (fabs(max - val) < _tolerance)
                return _slots_per_dim - 1;
            else if (fabs(min - val) < _tolerance)
                return 0;
            else {
                double pas = (fabs(max - min) / _slots_per_dim);
                int retorn = (int)((val - min) / pas);

                if (retorn >= _slots_per_dim)
                    return _slots_per_dim - 1;
                else if (retorn < 0)
                    return 0;
                else
                    return retorn;
            }
        }
    };

    ////////// GLOBALS (REMOVE)
    // Calculates numSimplification trajectory simplifications for one
    // trajectory
    // TODO: improve the epsilon to apply to more types of trajectories,
    // TODO: or track the simplification coarseness and adapt epsilon based on
    // it.
    void make_simplifications(ext_trajectory &t, double diagonal,
                              agarwal_simplification &agarwal, int size) {
        double targets[4] = {.07, .19, .24, .32};

        size_t targetCounts[4];
        for (int i = 0; i < 4; i++) {
            targetCounts[i] = static_cast<size_t>(t.size() * targets[i]);
            targetCounts[i] = std::max<size_t>(20, targetCounts[i]);
        }
        targetCounts[0] = std::min<size_t>(
            18, targetCounts[0]);  // start simple in case dihash is useless

        // double diag = t.bounding_box.get_diagonal();
        double lowerBound = diagonal / 100000;
        double upperBound = diagonal / 2;
        // int numIterations = 10;

        for (int i = 0; i < size; i++) {
            int tries = 0;
            double newUpperbound = 0.0;
            binary_double_search(
                [&](double value) -> int {
                    newUpperbound = value;
                    simplification simp = agarwal.simplify(t, value, _dist2);
                    tries++;
                    if (tries == 10) {  // TODO: test against numIterations???
                        t.add_simplification(std::move(simp));
                        return -1;
                    } else {
                        return simp.size() > targetCounts[i];
                    }
                },
                upperBound, lowerBound);
            upperBound = newUpperbound;
            // numIterations -= 2;
            // double ratio = simp->size/(double)t.size;
            _avgs_bb_ratio[i] += newUpperbound / diagonal;
        }
        _count++;

        /*
        for (int i = 0; i < size; i++) {
                simplification* ts = algo.agarwal.simplify(t,
        simplification_eps / (6 * (i + .25)));
                int diff = ts->size - t.simplifications[i]->size;
                avgs[i] += diff;
        }
        double avg = avgs[2] / (double)_count;
        std::cout << "avg " << avg << "\n";
        */

        // compile portals
        for (int i = 0; i < size; i++) {
            for (portal &p : t.simplifications[i].portals) {
                // check if it is a useful portal
                if (p.destination - p.source != 1) {
                    // check if it is not a duplicate
                    bool found = false;
                    for (portal &q : t.simp_portals[p.source]) {
                        if (q.destination == p.destination) {
                            found = true;
                        }
                    }
                    if (!found) {
                        t.simp_portals[p.source].push_back(p);
                    }
                }
            }
        }
        // sort portals from small to large
        for (auto &pair : t.simp_portals) {
            std::vector<portal> &k = pair.second;
            std::sort(k.begin(), k.end());
        }
    }

    // Calculates numSimplification trajectory simplifications for one
    // trajectory
    // TODO: improve the epsilon to apply to more types of trajectories,
    // TODO: or track the simplification coarseness and adapt epsilon based on
    // it.
    void make_source_simplifications(ext_trajectory &t, ext_trajectory &source,
                                     double diagonal,
                                     progressive_agarwal &agarwal_prog,
                                     int size) {
        for (int i = 0; i < size; i++) {
            double eps = diagonal * (_avgs_bb_ratio[i] / _count);
            t.add_simplification(agarwal_prog.simplify(t, source, eps, _dist2));
        }

        // compile portals
        for (int i = 0; i < size; i++) {
            for (portal &p : t.simplifications[i].portals) {
                // check if it is a useful portal
                if (p.destination - p.source != 1) {
                    // check if it is not a duplicate
                    bool found = false;
                    for (portal &q : t.simp_portals[p.source]) {
                        if (q.destination == p.destination) {
                            found = true;
                        }
                    }
                    if (!found) {
                        t.simp_portals[p.source].push_back(p);
                    }
                }
            }
        }

        // sort portals from small to large
        for (auto &pair : t.simp_portals) {
            std::vector<portal> &k = pair.second;
            std::sort(k.begin(), k.end());
        }
    }

    // pre-processing steps
    // --------------------------------------------------------------
    // Query step. Does rangequeries for start/endpoints of dataset. Adds all
    // found trajectories
    // to candidates.
    void collect_di_hash_points(
        const ext_trajectory &query_trajectory, double query_delta,
        squared_distance &dist2,
        const std::function<void(const ext_trajectory &)> &emit) {
        const auto &start = query_trajectory.front();
        const auto &end = query_trajectory.back();
        _di_hash.neighbors_with_callback(start, end, query_delta,
                                         _ext_trajectories, dist2, emit);
    }

    // Preprocessing step. Inserts start and endpoints in a regular grid so they
    // can be used
    // for range queries later.
    void add_pts_to_di_hash() {
        _di_hash.init(_bounding_box);
        for (ext_trajectory &t : _ext_trajectories) {
            if (t.size() > 1) {
                _di_hash.add_point(t.front(), t, true);
                _di_hash.add_point(t.back(), t, false);
            }
        }
    }

    void make_simplifications(ext_trajectory &t,
                              agarwal_simplification &agarwal) {
        make_simplifications(t, t.bounding_box.get_diagonal(), agarwal,
                             NUM_SIMPLIFICATIONS);
    }

    // Preprocessing step. Calculates simplifications for all trajectories in
    // the
    // dataset
    void construct_simplifications() {
#ifdef ENABLE_MULTITHREADING
        const size_t numWorkers =
            std::thread::hardware_concurrency();  // worker threads == number of
                                                  // logical cores
        if (numWorkers > 1 && _ext_trajectories.size() >= 100) {
            // use multiple threads
            std::vector<std::thread> simplificationThreads;
            std::vector<bounding_box_t> bboxes(numWorkers);

            // start threads
            const size_t trajsPerThread = _ext_trajectories.size() / numWorkers;
            for (size_t i = 0; i < numWorkers - 1; i++) {
                simplificationThreads.emplace_back(
                    &spatial_hash::simplification_worker, this,
                    i * trajsPerThread, (i + 1) * trajsPerThread, &(bboxes[i]));
            }
            simplificationThreads.emplace_back(
                &spatial_hash::simplification_worker, this,
                (numWorkers - 1) * trajsPerThread, _ext_trajectories.size(),
                &(bboxes[numWorkers - 1]));

            // join threads
            for (size_t i = 0; i < numWorkers; ++i) {
                simplificationThreads[i].join();
                _bounding_box.add_point(bboxes[i].minx, bboxes[i].miny);
                _bounding_box.add_point(bboxes[i].maxx, bboxes[i].maxy);
            }
            return;
        }
#endif

        // otherwise, use single thread
        simplification_worker(0, _ext_trajectories.size(), &_bounding_box);
    }

    void simplification_worker(size_t start_idx, size_t end_idx,
                               bounding_box_t *bbox) {
        agarwal_simplification agarwal;
        for (size_t i = start_idx; i < end_idx; ++i) {
            ext_trajectory &t = _ext_trajectories[i];
            if (t.size() > 1) {
                bbox->add_point(t.bounding_box.minx, t.bounding_box.miny);
                bbox->add_point(t.bounding_box.maxx, t.bounding_box.maxy);
                make_simplifications(t, agarwal);
            }
        }
    }

    /*
            PRUNING STRATEGIES
            ==================================
    */
    template <typename output_func>
    void return_trajectory(const ext_trajectory &et,
                           output_func &output) const {
        const trajectory &t = _src_trajectories[et.index];
        output(t);
    }

    // Query step. For each trajectory T in the dataset and query trajectory Q,
    // this step
    // compares successive simplifications of T and Q with continuous decision
    // frechet.
    // Each comparison can result in YES, NO, or MAYBE.
    // YES   -> remove from candidates, add to results
    // NO    -> remove from candidates
    // MAYBE -> try next simplification, if none are left, continue to next
    // algorithm step
    template <typename output_func>
    void prune_with_simplifications(
        cfdq_shortcuts &cdfqs, const ext_trajectory &query_trajectory,
        double query_delta, const ext_trajectory &t,
        const std::function<void(const ext_trajectory &)> &maybe,
        output_func &output) {
        bool broke = false;
        for (int i = 0; i < NUM_SIMPLIFICATIONS; ++i) {
            double decisionEpsilonLower =
                query_delta -
                query_trajectory.simplifications[i].simplification_eps -
                t.simplifications[i].simplification_eps;

            double decisionEpsilonUpper =
                query_delta +
                query_trajectory.simplifications[i].simplification_eps +
                t.simplifications[i].simplification_eps;

            double dist = equal_time_distance(
                t.simplifications[i], query_trajectory.simplifications[i]);

            if (dist < decisionEpsilonLower) {
                return_trajectory(t, output);
                broke = true;
                break;
            }

            if (decisionEpsilonLower > 0) {
                bool r = cdfqs.calculate(
                    query_trajectory.simplifications[i], t.simplifications[i],
                    decisionEpsilonLower, query_delta, _dist2);
                if (r) {
                    return_trajectory(t, output);
                    broke = true;
                    break;
                }
            }
            if (decisionEpsilonUpper > 0) {
                bool r = cdfqs.calculate(
                    query_trajectory.simplifications[i], t.simplifications[i],
                    decisionEpsilonUpper, query_delta, _dist2);
                if (!r) {
                    broke = true;
                    break;
                }
            }
        }
        if (!broke) {
            maybe(t);
        }
    }

    // Query step. Uses equal time distance as an upperbound for the actual
    // frechet distance
    // If ETD(P, Q) <= query_delta then CDF(P,Q) <= query_delta. With P in
    // dataset
    // and Q query trajectory.
    template <typename output_func>
    void prune_with_equal_time(
        const ext_trajectory &query_trajectory, double query_delta,
        const ext_trajectory &t,
        const std::function<void(const ext_trajectory &)> &maybe,
        output_func &output) {
        double dist = equal_time_distance(t, query_trajectory);
        if (dist < query_delta) {
            return_trajectory(t, output);
        } else {
            maybe(t);
        }
    }

    // Query step. The final step for each query is to do a full decision
    // frechet
    // computation.
    // This step contains no additional smart optimization, and so is very slow.
    template <typename output_func>
    void prune_with_decision_frechet(cfdq_shortcuts &cdfqs,
                                     const ext_trajectory &query_trajectory,
                                     double query_delta,
                                     const ext_trajectory &t,
                                     squared_distance &dist2,
                                     output_func &output) {
        if (cdfqs.calculate(query_trajectory, t, query_delta, dist2)) {
            return_trajectory(t, output);
        }
    }

    // -----------------------------------------

    // Only computes reachable part of freespace diagram, uses shortcuts to skip
    // columns in reachable part
    class cfdq_shortcuts {
        struct q_entry {
            int start_row_index;
            int end_row_index;
            double lowest_right;
        };

        std::vector<q_entry> queue[2];
        int queue_size[2];

        double compute_segment_frechet(const portal &p, int q,
                                       const vertices &p_array,
                                       const vertices &q_array,
                                       squared_distance &dist2) {
            const point &pstart = p_array[p.source];
            const point &pend = p_array[p.destination];
            const point &qstart = q_array[q];
            const point &qend = q_array[q];
            return std::sqrt(
                std::max(dist2(pstart, qstart), dist2(pend, qend)));
        }

       public:
        int numRows = 0;
        bool calculate(const vertices &P, const vertices &Q, size_t offset_p,
                       size_t offset_q, size_t size_p, size_t size_q,
                       double query_delta, double base_query_delta,
                       const std::map<int, std::vector<portal>> &portals,
                       squared_distance &dist2) {
            double startDist = dist2(P[offset_p], Q[offset_q]);
            double endDist = dist2(P[size_p - 1], Q[size_q - 1]);
            if (startDist > query_delta * query_delta ||
                endDist > query_delta * query_delta)
                return false;
            if (size_p <= offset_p + 1 || size_q <= offset_q + 1)
                return false;  // TODO: do we need this? // WM: added offset

            int first = 0;
            int second = 1;

            range Rf;
            range Tf;
            portal choice;

            // ensure queue capacity
            // WM: added offsets
            size_t max = size_p - offset_p;
            if (size_q - offset_q > max) max = size_q - offset_q;
            if (queue[0].size() < max) {
                for (size_t i = queue[0].size(); i < max; i++) {
                    queue[0].push_back({0, 0});
                    queue[1].push_back({0, 0});
                }
            }

            // setup
            // WM: this is always free space! by check on startDist
            // bool LFree = compute_interval(Q[0], P[0], P[1], query_delta, Rf);
            queue[first][0].start_row_index = 0;
            queue[first][0].end_row_index = 0;
            queue[first][0].lowest_right = 0;

            queue_size[first] = 1;
            queue_size[second] = 0;

            // For each column
            for (int column = offset_q; column < static_cast<int>(size_q) - 1;
                 column++) {
                if (queue_size[first] == 0) {
                    // nothing reachable anymore
                    return false;
                }
                queue_size[second] = 0;
                int row = queue[first][0].start_row_index;
                int qIndex = 0;
                // while there's reachable cells left in the queue
                while (qIndex < queue_size[first]) {
                    double left_most_top = 2;
                    // start at reachable cell at the head of the queue, and
                    // continue
                    // until
                    // reachability cannot propagate, consuming the queue as we
                    // progress
                    do {
                        // tracks whether we overshoot the queue
                        bool outsideQueue = qIndex >= queue_size[first];
                        // Right edge stored in Rf, RFree = false means not free
                        bool RFree = compute_interval(
                            Q[column + 1], P[row], P[row + 1], query_delta, Rf);
                        if (RFree) {
                            if (left_most_top <= 1) {
                                double newLR = Rf.start;
                                if (Rf.is_complete() &&
                                    queue_size[second] > 0 &&
                                    queue[second][queue_size[second] - 1]
                                            .end_row_index == row - 1) {
                                    // complete reachable right means increase
                                    // previous queue
                                    // entry to span this cell
                                    queue[second][queue_size[second] - 1]
                                        .end_row_index = row;
                                } else {
                                    // push to queue
                                    queue[second][queue_size[second]]
                                        .start_row_index = row;
                                    queue[second][queue_size[second]]
                                        .end_row_index = row;
                                    queue[second][queue_size[second]]
                                        .lowest_right = newLR;
                                    queue_size[second]++;
                                }
                            } else {
                                // WM: think you should be checking row here as
                                // well
                                if (!outsideQueue &&
                                    row >=
                                        queue[first][qIndex].start_row_index &&
                                    row <= queue[first][qIndex].end_row_index) {
                                    if (!(row ==
                                              queue[first][qIndex]
                                                  .start_row_index &&
                                          queue[first][qIndex].lowest_right >
                                              Rf.end)) {
                                        double prevR =
                                            row ==
                                                    queue[first][qIndex]
                                                        .start_row_index
                                                ? queue[first][qIndex]
                                                      .lowest_right
                                                : 0.0;
                                        double newLR =
                                            std::max(prevR, Rf.start);
                                        if (Rf.is_complete() && newLR == 0.0 &&
                                            queue_size[second] > 0 &&
                                            queue[second]
                                                 [queue_size[second] - 1]
                                                     .end_row_index ==
                                                row - 1) {
                                            // complete reachable right means
                                            // increase previous queue
                                            // entry to span this cell
                                            queue[second]
                                                 [queue_size[second] - 1]
                                                     .end_row_index = row;
                                        } else {
                                            // push to queue
                                            queue[second][queue_size[second]]
                                                .start_row_index = row;
                                            queue[second][queue_size[second]]
                                                .end_row_index = row;
                                            queue[second][queue_size[second]]
                                                .lowest_right = newLR;
                                            queue_size[second]++;
                                        }
                                    }
                                }
                            }
                        }
                        // Top edge stored in Tf, TFree = false means not free
                        bool TFree =
                            compute_interval(P[row + 1], Q[column],
                                             Q[column + 1], query_delta, Tf);
                        if (!outsideQueue &&
                            row <= queue[first][qIndex].end_row_index &&
                            row >= queue[first][qIndex].start_row_index) {
                            if (row == queue[first][qIndex].end_row_index) {
                                // consume the first queue
                                qIndex++;
                            }
                            if (TFree) {
                                left_most_top = Tf.start;
                            } else {
                                left_most_top = 2;
                            }
                        } else if (TFree && left_most_top <= Tf.end) {
                            left_most_top = std::max(left_most_top, Tf.start);
                        } else {
                            left_most_top = 2;
                        }
                        // try and jump
                        if (!outsideQueue && queue_size[second] > 0 &&
                            queue[second][queue_size[second] - 1]
                                    .end_row_index == row &&
                            Rf.end == 1.0) {
                            // jump-off point possible
                            // check if minimum jump distance is big enough
                            int gapSize = queue[first][qIndex].end_row_index -
                                          queue[first][qIndex].start_row_index;
                            if (gapSize > 1) {
                                const auto iter = portals.find(row);
                                if (iter != portals.cend()) {
                                    const std::vector<portal> &ports =
                                        iter->second;
                                    choice.source = -1;
                                    for (const portal &p : ports) {
                                        // int jumpSize = p.destination -
                                        // p.source;
                                        // check if jump within range
                                        if (p.destination <=
                                            queue[first][qIndex]
                                                .end_row_index) {
                                            // check if jump distance fits
                                            double segmentFrechet =
                                                compute_segment_frechet(
                                                    p, column, P, Q, dist2);
                                            if (segmentFrechet + p.distance <=
                                                base_query_delta) {
                                                choice = p;
                                            }
                                        } else {
                                            break;
                                        }
                                    }
                                    // JUMP!
                                    if (choice.source != -1) {
                                        row = choice.destination -
                                              1;  // - 1 to counter ++ later
                                        queue[second][queue_size[second] - 1]
                                            .end_row_index = row;
                                    }
                                }
                            }
                        }
                        // propagated reachability by one cell, so look at next
                        // row
                        row++;
                        numRows++;
                    } while (left_most_top <= 1 &&
                             row < static_cast<int>(size_p) - 1);
                }

                // swap first and second column
                int temp = first;
                first = second;
                second = temp;
            }

            int endIndex = queue_size[first] - 1;
            if (endIndex == -1) return false;
            bool exit = queue[first][endIndex].start_row_index ==
                            static_cast<int>(size_p) - 2 &&
                        queue[first][endIndex].lowest_right <= 1;
            return exit || (queue[first][endIndex].end_row_index ==
                                static_cast<int>(size_p) - 2 &&
                            queue[first][endIndex].start_row_index !=
                                static_cast<int>(size_p) - 2);
        }

        bool calculate(const ext_trajectory &P, const ext_trajectory &Q,
                       double query_delta, double base_query_delta,
                       squared_distance &dist2) {
            return calculate(P.points, Q.points, 0, 0, P.size(), Q.size(),
                             query_delta, base_query_delta, P.simp_portals,
                             dist2);
        }

        bool calculate(const ext_trajectory &P, const ext_trajectory &Q,
                       double query_delta, squared_distance &dist2) {
            return calculate(P.points, Q.points, 0, 0, P.size(), Q.size(),
                             query_delta, query_delta, P.simp_portals, dist2);
        }
    };

    // Does binary search on integer range (lowerbound, upperbound),
    // accepts lambda function returning whether the given search index
    // satisfies the search criterion.

    void binary_double_search(const std::function<int(double)> &f,
                              double upperbound, double lowerbound) {
        double rangeLength = upperbound - lowerbound;
        double avg = lowerbound + (rangeLength) / 2;
        int result = f(avg);
        if (result == 1) {
            binary_double_search(f, upperbound, avg);
        } else if (result == 0) {
            binary_double_search(f, avg, lowerbound);
        } else {
            return;
        }
    }

    // Does double & search on integer range (lowerbound, upperbound),
    // accepts lambda function returning whether the given search index
    // satisfies the search criterion.
    static int double_n_search(const std::function<bool(int)> &f, int start,
                               int end, int double_n_search_base,
                               double exponent_step) {
        int k = start;
        int prevk = start;
        int iteration = 0;
        while (true) {
            // double
            if (k > end - 1) {
                k = end - 1;
            }
            bool epsValid = f(k);
            if (!epsValid) {
                // binary search
                k = binary_int_search(f, k, prevk);
                return k;
            } else {
                if (k == end - 1) {
                    return k;
                }
                prevk = k;
                k += (int)floor(
                    pow(double_n_search_base, exponent_step * iteration));
                iteration++;
            }
        }
    }

    static int binary_int_search(const std::function<bool(int)> &f,
                                 int upperbound, int lowerbound) {
        int rangeLength = upperbound - lowerbound;
        if (rangeLength <= 1) {
            return lowerbound;
        }
        int middle = lowerbound + (rangeLength) / 2;
        bool result = f(middle);
        if (result) {
            return binary_int_search(f, upperbound, middle);
        } else {
            return binary_int_search(f, middle, lowerbound);
        }
    }

    // This class contains all logic need to compute a simplification
    // of any input trajectory using agarwal simplification with
    // double & search. The algorithm uses EqualTimeDistance.h to
    // satisfy the agarwal constraints.
    class agarwal_simplification {
       public:
        simplification simplify(ext_trajectory &t, double simplification_eps,
                                squared_distance &dist2) {
            simplification simplified(simplification_eps);

            const vertices &P = t.points;
            simplified.points.push_back(P[0]);
            simplified.distances.push_back(0);
            simplified.totals.push_back(0);

            int simp_size = 1;

            int rangeStart = 1;
            int prevk = 0;
            while (true) {
                int k = find_last_frechet_match(
                    P, simplified.points, t.totals, simplified.totals,
                    t.distances, simplified.distances,
                    static_cast<int>(simplified.points.size()), rangeStart,
                    static_cast<int>(t.size()), prevk, simplification_eps,
                    dist2);
                simp_size++;
                simplified.points[simp_size - 1] = P[k];
                if (k == static_cast<int>(t.size()) - 1) {
                    break;
                }
                prevk = k;
                rangeStart = k + 1;
            }

            return simplified;
        }

       private:
        // Finds index k of last vertex v that still satisfies
        // the simplification epsilon
        int find_last_frechet_match(
            const vertices &P, vertices &simp, std::vector<double> &ptotals,
            std::vector<double> &simptotals, std::vector<double> &pdists,
            std::vector<double> &simpdists, int simp_size, int start, int end,
            int prevk, double epsilon, squared_distance &dist2) {
            simp.push_back(P[0]);
            simpdists.push_back(0);
            simptotals.push_back(0);
            // Use lambda's to easily double & search the function from (start)
            // to
            // (end)
            constexpr int double_n_search_base = 2;
            constexpr double exponent_step = 1.0;
            return double_n_search(
                [&](int index) -> bool {
                    simp[simp_size] = P[index];
                    simpdists[simp_size] =
                        std::sqrt(dist2(simp[simp_size], simp[simp_size - 1]));
                    simptotals[simp_size] =
                        simptotals[simp_size - 1] + simpdists[simp_size];
                    double dist = equal_time_distance(
                        P, simp, ptotals, simptotals, pdists, simpdists,
                        index + 1, simp_size + 1, prevk, simp_size - 1);
                    return dist <= epsilon;
                },
                start, end, double_n_search_base, exponent_step);
        }
    };

    // This class contains all logic need to compute a simplification
    // of any input trajectory using agarwal simplification with
    // double & search. The algorithm uses EqualTimeDistance.h to
    // satisfy the agarwal constraints.
    class progressive_agarwal {
       public:
        simplification simplify(ext_trajectory &parent,
                                ext_trajectory &source_trajectory,
                                double simplification_eps,
                                squared_distance &dist2) {
            simplification simplified(simplification_eps);

            const vertices &P = parent.points;
            simplified.points.push_back(P[0]);
            simplified.distances.push_back(0);
            simplified.totals.push_back(0);
            simplified.sourceIndex.push_back(0);

            int simp_size = 1;

            int rangeStart = 1;
            int prevk = 0;
            while (true) {
                int k = find_last_frechet_match(
                    P, simplified.points, parent.totals, simplified.totals,
                    parent.distances, simplified.distances,
                    static_cast<int>(simplified.points.size()), rangeStart,
                    static_cast<int>(parent.size()), prevk, simplification_eps,
                    dist2, parent.sourceIndex, source_trajectory,
                    simplified.portals);
                simp_size++;
                simplified.points[simp_size - 1] = P[k];
                simplified.sourceIndex.push_back(parent.sourceIndex[k]);
                if (k == static_cast<int>(parent.size()) - 1) {
                    break;
                }
                prevk = k;
                rangeStart = k + 1;
            }

            return simplified;
        }

       private:
        // Finds index k of last vertex v that still satisfies
        // the simplification epsilon
        int find_last_frechet_match(
            const vertices &P, vertices &simp, std::vector<double> &ptotals,
            std::vector<double> &simptotals, std::vector<double> &pdists,
            std::vector<double> &simpdists, int simp_size, int start, int end,
            int prevk, double epsilon, squared_distance &dist2,
            std::vector<size_t> &parent_source_indices,
            ext_trajectory &source_trajectory, std::vector<portal> &portals) {
            simp.push_back(P[0]);
            simpdists.push_back(0);
            simptotals.push_back(0);
            // Use lambda's to easily double & search the function from (start)
            // to
            // (end)
            constexpr int double_n_search_base = 2;
            constexpr double exponent_step = 1.0;
            return double_n_search(
                [&](int index) -> bool {
                    simp[simp_size] = P[index];
                    simpdists[simp_size] =
                        std::sqrt(dist2(simp[simp_size], simp[simp_size - 1]));
                    simptotals[simp_size] =
                        simptotals[simp_size - 1] + simpdists[simp_size];
                    size_t end = 0;
                    size_t start = parent_source_indices[prevk];
                    if (static_cast<size_t>(index + 1) >=
                        parent_source_indices.size()) {
                        end = source_trajectory.size();
                    } else {
                        end = parent_source_indices[index + 1];
                    }
                    double dist = equal_time_distance(
                        source_trajectory.points, simp,
                        source_trajectory.totals, simptotals,
                        source_trajectory.distances, simpdists,
                        static_cast<int>(end), simp_size + 1,
                        static_cast<int>(start), simp_size - 1);
                    portal p;
                    p.source = prevk;
                    p.destination = index;
                    p.distance = dist;
                    portals.push_back(p);
                    return dist <= epsilon;
                },
                start, end, double_n_search_base, exponent_step);
        }
    };
};

}  // namespace bddm
}  // namespace detail
}  // namespace frechetrange

#endif  // TUE_INC
