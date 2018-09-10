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
    typedef std::vector<point> Vertices;
    struct BoundingBox;
    struct Portal;
    struct ExtTrajectory;
    struct TrajectorySimplification;

    class DiHash;
    class AgarwalSimplification;
    class ProgressiveAgarwal;
    class CDFQShortcuts;

   public:
    spatial_hash(const squared_distance &dist2 = squared_distance())
        : _dist2(dist2),
          _srcTrajectories(),
          _extTrajectories(),
          _boundingBox(),
          _diHash(SLOTS_PER_DIMENSION, TOLERANCE),
          _count(0),
          _avgsBBRatio{0.0, 0.0, 0.0, 0.0} {}

    size_t size() const { return _srcTrajectories.size(); }

    void insert(const trajectory &t) {
        if (t.size() > 0) {
            _srcTrajectories.push_back(t);
            _extTrajectories.emplace_back(t, _extTrajectories.size(), _dist2);
        }
    }

    void insert(trajectory &&t) {
        if (t.size() > 0) {
            _extTrajectories.emplace_back(t, _extTrajectories.size(), _dist2);
            _srcTrajectories.emplace_back(std::move(t));
        }
    }

    // Does all needed preprocessing for the given dataset
    void build_index() {
        constructSimplifications();
        addPtsToDiHash();
    }

    std::vector<const trajectory *> range_query(const trajectory &query,
                                                double distanceThreshold) {
        std::vector<const trajectory *> resultSet;
        auto pushBackResult = [&resultSet](const trajectory &t) {
            resultSet.push_back(&t);
        };
        range_query(query, distanceThreshold, pushBackResult);
        return resultSet;
    }

    /**
    * @param output Supports the method operator()(const trajectory &)
    *               to output the result trajectories.
    */
    template <typename OutputFunctional>
    void range_query(const trajectory &query, double distanceThreshold,
                     OutputFunctional &output) {
        ExtTrajectory queryTrajectory(query, 0, _dist2);

        ProgressiveAgarwal agarwalProg;
        double diagonal = queryTrajectory.boundingBox.getDiagonal();
        makeSourceSimplificationsForTrajectory(queryTrajectory, queryTrajectory,
                                               diagonal, agarwalProg,
                                               NUM_SIMPLIFICATIONS);
        for (int i = 1; i < NUM_SIMPLIFICATIONS; i++) {
            makeSourceSimplificationsForTrajectory(
                queryTrajectory.simplifications[i], queryTrajectory, diagonal,
                agarwalProg, i - 1);
        }

        CDFQShortcuts cdfqs;
        collectDiHashPoints(
            queryTrajectory, distanceThreshold, _dist2,
            [&](const ExtTrajectory &t) {
                pruneWithSimplifications(
                    cdfqs, queryTrajectory, distanceThreshold, t,
                    [&](const ExtTrajectory &t) {
                        pruneWithEqualTime(
                            queryTrajectory, distanceThreshold, t,
                            [&](const ExtTrajectory &t) {
                                pruneWithDecisionFrechet(cdfqs, queryTrajectory,
                                                         distanceThreshold, t,
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

    std::vector<trajectory> _srcTrajectories;
    std::vector<ExtTrajectory> _extTrajectories;
    BoundingBox _boundingBox;
    DiHash _diHash;

    int _count;
    double _avgsBBRatio[4];

    /* DATA STRUCTURES */

    // Represents a boundingbox around a trajectory or group of trajectories
    struct BoundingBox {
        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::max();

        double maxx = std::numeric_limits<double>::min();
        double maxy = std::numeric_limits<double>::min();

        void addPoint(const point &p) {
            addPoint(get_coordinate::template get<0>(p),
                     get_coordinate::template get<1>(p));
        }

        void addPoint(double x, double y) {
            if (x < minx) minx = x;
            if (y < miny) miny = y;

            if (x > maxx) maxx = x;
            if (y > maxy) maxy = y;
        }

        double getDiagonal() {
            double width = maxx - minx;
            double height = maxy - miny;
            return std::sqrt(width * width + height * height);
        }
    };

    struct Portal {
        int source;
        int destination;
        double distance;

        bool operator<(const Portal &p) const {
            return destination < p.destination;
        };
    };

    struct ExtTrajectory {
        size_t index;
        Vertices vertices;
        std::vector<double> distances;  // distance between vertices
        std::vector<double> totals;     // total length at vertex x
        std::vector<size_t>
            sourceIndex;  // mapping back to the source trajectory if
                          // this is a simplification
        BoundingBox boundingBox;

        std::map<int, std::vector<Portal>> simpPortals;  // freespace jumps
        std::vector<TrajectorySimplification> simplifications;

        ExtTrajectory() = default;
        ExtTrajectory(ExtTrajectory &&) = default;

        ExtTrajectory(const trajectory &t, size_t idx, squared_distance &dist2)
            : index(idx) {
            vertices.reserve(t.size());
            vertices.push_back(t[0]);
            distances.push_back(0.0);
            totals.push_back(0.0);
            sourceIndex.push_back(0);
            boundingBox.addPoint(t[0]);

            for (size_t i = 1; i < t.size(); ++i) {
                const point &p = t[i];
                double x = get_coordinate::template get<0>(p);
                double y = get_coordinate::template get<1>(p);

                // ignore NaN and duplicates
                const point &prev = vertices.back();
                if (std::isfinite(x) && std::isfinite(y) &&
                    (x != get_coordinate::template get<0>(prev) ||
                     y != get_coordinate::template get<1>(prev))) {
                    vertices.push_back(p);
                    distances.push_back(
                        std::sqrt(dist2(prev, vertices.back())));
                    totals.push_back(totals.back() + distances.back());
                    sourceIndex.push_back(sourceIndex.size());
                    boundingBox.addPoint(p);
                }
            }
        }

        size_t size() const { return vertices.size(); }

        const point &front() const { return vertices[0]; }

        const point &back() const { return vertices.back(); }

        void addSimplification(TrajectorySimplification &&simp) {
            simplifications.emplace_back(std::move(simp));
        }
    };

    struct TrajectorySimplification : public ExtTrajectory {
        double simplificationEpsilon;
        std::vector<Portal> portals;

        TrajectorySimplification() = default;
        TrajectorySimplification(TrajectorySimplification &&) = default;
        TrajectorySimplification(double eps)
            : simplificationEpsilon(eps), portals() {}
    };

    static double equalTimeDistance(const ExtTrajectory &p,
                                    const ExtTrajectory &q) {
        return equalTimeDistance(p.vertices, q.vertices, p.totals, q.totals,
                                 p.distances, q.distances, p.size(), q.size(),
                                 0, 0);
    }

    // Implements Equal Time Distance algorithm between two trajectories.
    // The ETD algorithm computes an approximation of frechet distance by
    // taking the 'dog leash' length when traversing two trajectories at
    // the same speed. Used by agarwal and simplification step.
    // If found relevant, this function can be optimized using SIMD instructions
    static double equalTimeDistance(const Vertices &pverts,
                                    const Vertices &qverts,
                                    const std::vector<double> &ptotals,
                                    const std::vector<double> &qtotals,
                                    const std::vector<double> &pdistances,
                                    const std::vector<double> &qdistances,
                                    size_t psize, size_t qsize, size_t pstart,
                                    size_t qstart) {
        double pdistOffset = ptotals[pstart];
        double qdistOffset = qtotals[qstart];
        double pdist = ptotals[psize - 1] - pdistOffset;
        double qdist = qtotals[qsize - 1] - qdistOffset;
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
            double nextDistP = ptotals[p_ptr] - pdistOffset - posP;
            double nextDistQ = qtotals[q_ptr] - qdistOffset - posQ;

            if (p_ptr == psize - 1)
                nextDistP = std::numeric_limits<double>::max();
            if (q_ptr == qsize - 1)
                nextDistQ = std::numeric_limits<double>::max();

            if (nextDistP * pscale < nextDistQ) {  // treat P first
                p_ptx = get_coordinate::template get<0>(pverts[p_ptr]);
                p_pty = get_coordinate::template get<1>(pverts[p_ptr]);
                position = (ptotals[p_ptr] - pdistOffset) / pdist;
                double scale =
                    (position * qdist - (qtotals[q_ptr - 1] - qdistOffset)) /
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

                position = (qtotals[q_ptr] - qdistOffset) / qdist;
                double scale =
                    (position * pdist - (ptotals[p_ptr - 1] - pdistOffset)) /
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

    struct Range {
        double start = 0.0;
        double end = 0.0;

        bool isComplete() const { return start == 0.0 && end == 1.0; }
    };

    static bool computeInterval(const point &a, const point &b1,
                                const point &b2, double eps, Range &r) {
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
    class DiHash {
        struct Vertex {
            const point &p;
            const ExtTrajectory &t;
            bool isStart;
        };

        const int
            _slotsPerDimension;  // number of equally spaced subdivisions in
                                 // each dimension
        std::vector<std::vector<Vertex>> _elements;
        double
            _tolerance;  // tolerance to prevent numerical representation errors
        double
            _limits[2]
                   [2];  // two rows (x,y) and two columns(min,max), keeps the
                         // information on limits in each dimension, for 3D,
                         // just add one row.

       public:
        DiHash(int numC, double tol)
            : _slotsPerDimension(numC),
              _elements(numC * numC),
              _tolerance(tol) {}

        void init(const BoundingBox &boundingBox) {
            // First, locate numerichal limits
            // initialize limits at a three-component two doubles (max and min)
            // vector
            // set initial values for extreme values
            _limits[0][0] = boundingBox.minx;  // min in X
            _limits[0][1] = boundingBox.maxx;  // max in X
            _limits[1][0] = boundingBox.miny;  // min in Y
            _limits[1][1] = boundingBox.maxy;  // max in Y
        }

        void addPoint(const point &p, const ExtTrajectory &t, bool isStart) {
            int xSlot =
                findSlot(get_coordinate::template get<0>(p), 'x', false);
            int ySlot =
                findSlot(get_coordinate::template get<1>(p), 'y', false);
            _elements[xSlot * _slotsPerDimension + ySlot].push_back(
                Vertex{p, t, isStart});
        }

        // Return neighbors at range strictly less than distance for a given
        // point,
        // does not return the query point if present
        void neighborsWithCallback(
            const point &start, const point &end, double eps,
            std::vector<ExtTrajectory> &trajectories, squared_distance &dist2,
            const std::function<void(const ExtTrajectory &)> &emit) {
            auto limitsX =
                slotsTouched(get_coordinate::template get<0>(start), eps, 'x');
            auto limitsY =
                slotsTouched(get_coordinate::template get<1>(start), eps, 'y');
            for (int i = limitsX.first; i <= limitsX.second; i++) {
                for (int j = limitsY.first; j <= limitsY.second; j++) {
                    for (const Vertex &pActual :
                         _elements[i * _slotsPerDimension + j]) {
                        if (pActual.isStart) {
                            double distSQ = dist2(start, pActual.p);
                            if (distSQ < eps * eps) {
                                const ExtTrajectory &t = pActual.t;
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
        std::pair<int, int> slotsTouched(double center, double eps, char type) {
            double min = center - eps;
            double max = center + eps;
            // first position of the following vector is for the minimum slot
            // affected
            // and second for the maximum
            return std::make_pair(findSlot(min, type, true),
                                  findSlot(max, type, true));
        }

        int findSlot(double val, char type, bool allowOverflow) {
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
                return _slotsPerDimension - 1;
            else if (fabs(min - val) < _tolerance)
                return 0;
            else {
                double pas = (fabs(max - min) / _slotsPerDimension);
                int retorn = (int)((val - min) / pas);

                if (retorn >= _slotsPerDimension)
                    return _slotsPerDimension - 1;
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
    void makeSimplificationsForTrajectory(ExtTrajectory &t, double diagonal,
                                          AgarwalSimplification &agarwal,
                                          int size) {
        double targets[4] = {.07, .19, .24, .32};

        size_t targetCounts[4];
        for (int i = 0; i < 4; i++) {
            targetCounts[i] = static_cast<size_t>(t.size() * targets[i]);
            targetCounts[i] = std::max<size_t>(20, targetCounts[i]);
        }
        targetCounts[0] = std::min<size_t>(
            18, targetCounts[0]);  // start simple in case dihash is useless

        // double diag = t.boundingBox.getDiagonal();
        double lowerBound = diagonal / 100000;
        double upperBound = diagonal / 2;
        // int numIterations = 10;

        for (int i = 0; i < size; i++) {
            int tries = 0;
            double newUpperbound = 0.0;
            binaryDoubleSearch(
                [&](double value) -> int {
                    newUpperbound = value;
                    TrajectorySimplification simp =
                        agarwal.simplify(t, value, _dist2);
                    tries++;
                    if (tries == 10) {  // TODO: test against numIterations???
                        t.addSimplification(std::move(simp));
                        return -1;
                    } else {
                        return simp.size() > targetCounts[i];
                    }
                },
                upperBound, lowerBound);
            upperBound = newUpperbound;
            // numIterations -= 2;
            // double ratio = simp->size/(double)t.size;
            _avgsBBRatio[i] += newUpperbound / diagonal;
        }
        _count++;

        /*
        for (int i = 0; i < size; i++) {
                TrajectorySimplification* ts = algo.agarwal.simplify(t,
        simplificationEpsilon / (6 * (i + .25)));
                int diff = ts->size - t.simplifications[i]->size;
                avgs[i] += diff;
        }
        double avg = avgs[2] / (double)_count;
        std::cout << "avg " << avg << "\n";
        */

        // compile portals
        for (int i = 0; i < size; i++) {
            for (Portal &p : t.simplifications[i].portals) {
                // check if it is a useful portal
                if (p.destination - p.source != 1) {
                    // check if it is not a duplicate
                    bool found = false;
                    for (Portal &q : t.simpPortals[p.source]) {
                        if (q.destination == p.destination) {
                            found = true;
                        }
                    }
                    if (!found) {
                        t.simpPortals[p.source].push_back(p);
                    }
                }
            }
        }
        // sort portals from small to large
        for (auto &pair : t.simpPortals) {
            std::vector<Portal> &k = pair.second;
            std::sort(k.begin(), k.end());
        }
    }

    // Calculates numSimplification trajectory simplifications for one
    // trajectory
    // TODO: improve the epsilon to apply to more types of trajectories,
    // TODO: or track the simplification coarseness and adapt epsilon based on
    // it.
    void makeSourceSimplificationsForTrajectory(ExtTrajectory &t,
                                                ExtTrajectory &source,
                                                double diagonal,
                                                ProgressiveAgarwal &agarwalProg,
                                                int size) {
        for (int i = 0; i < size; i++) {
            double eps = diagonal * (_avgsBBRatio[i] / _count);
            t.addSimplification(agarwalProg.simplify(t, source, eps, _dist2));
        }

        // compile portals
        for (int i = 0; i < size; i++) {
            for (Portal &p : t.simplifications[i].portals) {
                // check if it is a useful portal
                if (p.destination - p.source != 1) {
                    // check if it is not a duplicate
                    bool found = false;
                    for (Portal &q : t.simpPortals[p.source]) {
                        if (q.destination == p.destination) {
                            found = true;
                        }
                    }
                    if (!found) {
                        t.simpPortals[p.source].push_back(p);
                    }
                }
            }
        }

        // sort portals from small to large
        for (auto &pair : t.simpPortals) {
            std::vector<Portal> &k = pair.second;
            std::sort(k.begin(), k.end());
        }
    }

    // pre-processing steps
    // --------------------------------------------------------------
    // Query step. Does rangequeries for start/endpoints of dataset. Adds all
    // found trajectories
    // to candidates.
    void collectDiHashPoints(
        const ExtTrajectory &queryTrajectory, double queryDelta,
        squared_distance &dist2,
        const std::function<void(const ExtTrajectory &)> &emit) {
        const auto &start = queryTrajectory.front();
        const auto &end = queryTrajectory.back();
        _diHash.neighborsWithCallback(start, end, queryDelta, _extTrajectories,
                                      dist2, emit);
    }

    // Preprocessing step. Inserts start and endpoints in a regular grid so they
    // can be used
    // for range queries later.
    void addPtsToDiHash() {
        _diHash.init(_boundingBox);
        for (ExtTrajectory &t : _extTrajectories) {
            if (t.size() > 1) {
                _diHash.addPoint(t.front(), t, true);
                _diHash.addPoint(t.back(), t, false);
            }
        }
    }

    void makeSimplificationsForTrajectory(ExtTrajectory &t,
                                          AgarwalSimplification &agarwal) {
        makeSimplificationsForTrajectory(t, t.boundingBox.getDiagonal(),
                                         agarwal, NUM_SIMPLIFICATIONS);
    }

    // Preprocessing step. Calculates simplifications for all trajectories in
    // the
    // dataset
    void constructSimplifications() {
#ifdef ENABLE_MULTITHREADING
        const size_t numWorkers =
            std::thread::hardware_concurrency();  // worker threads == number of
                                                  // logical cores
        if (numWorkers > 1 && _extTrajectories.size() >= 100) {
            // use multiple threads
            std::vector<std::thread> simplificationThreads;
            std::vector<BoundingBox> bboxes(numWorkers);

            // start threads
            const size_t trajsPerThread = _extTrajectories.size() / numWorkers;
            for (size_t i = 0; i < numWorkers - 1; i++) {
                simplificationThreads.emplace_back(
                    &spatial_hash::simplificationWorker, this,
                    i * trajsPerThread, (i + 1) * trajsPerThread, &(bboxes[i]));
            }
            simplificationThreads.emplace_back(
                &spatial_hash::simplificationWorker, this,
                (numWorkers - 1) * trajsPerThread, _extTrajectories.size(),
                &(bboxes[numWorkers - 1]));

            // join threads
            for (size_t i = 0; i < numWorkers; ++i) {
                simplificationThreads[i].join();
                _boundingBox.addPoint(bboxes[i].minx, bboxes[i].miny);
                _boundingBox.addPoint(bboxes[i].maxx, bboxes[i].maxy);
            }
            return;
        }
#endif

        // otherwise, use single thread
        simplificationWorker(0, _extTrajectories.size(), &_boundingBox);
    }

    void simplificationWorker(size_t startIdx, size_t endIdx,
                              BoundingBox *bbox) {
        AgarwalSimplification agarwal;
        for (size_t i = startIdx; i < endIdx; ++i) {
            ExtTrajectory &t = _extTrajectories[i];
            if (t.size() > 1) {
                bbox->addPoint(t.boundingBox.minx, t.boundingBox.miny);
                bbox->addPoint(t.boundingBox.maxx, t.boundingBox.maxy);
                makeSimplificationsForTrajectory(t, agarwal);
            }
        }
    }

    /*
            PRUNING STRATEGIES
            ==================================
    */
    template <typename OutputFunctional>
    void returnTrajectory(const ExtTrajectory &et,
                          OutputFunctional &output) const {
        const trajectory &t = _srcTrajectories[et.index];
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
    template <typename OutputFunctional>
    void pruneWithSimplifications(
        CDFQShortcuts &cdfqs, const ExtTrajectory &queryTrajectory,
        double queryDelta, const ExtTrajectory &t,
        const std::function<void(const ExtTrajectory &)> &maybe,
        OutputFunctional &output) {
        bool broke = false;
        for (int i = 0; i < NUM_SIMPLIFICATIONS; ++i) {
            double decisionEpsilonLower =
                queryDelta -
                queryTrajectory.simplifications[i].simplificationEpsilon -
                t.simplifications[i].simplificationEpsilon;

            double decisionEpsilonUpper =
                queryDelta +
                queryTrajectory.simplifications[i].simplificationEpsilon +
                t.simplifications[i].simplificationEpsilon;

            double dist = equalTimeDistance(t.simplifications[i],
                                            queryTrajectory.simplifications[i]);

            if (dist < decisionEpsilonLower) {
                returnTrajectory(t, output);
                broke = true;
                break;
            }

            if (decisionEpsilonLower > 0) {
                bool r = cdfqs.calculate(
                    queryTrajectory.simplifications[i], t.simplifications[i],
                    decisionEpsilonLower, queryDelta, _dist2);
                if (r) {
                    returnTrajectory(t, output);
                    broke = true;
                    break;
                }
            }
            if (decisionEpsilonUpper > 0) {
                bool r = cdfqs.calculate(
                    queryTrajectory.simplifications[i], t.simplifications[i],
                    decisionEpsilonUpper, queryDelta, _dist2);
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
    // If ETD(P, Q) <= queryDelta then CDF(P,Q) <= queryDelta. With P in dataset
    // and Q query trajectory.
    template <typename OutputFunctional>
    void pruneWithEqualTime(
        const ExtTrajectory &queryTrajectory, double queryDelta,
        const ExtTrajectory &t,
        const std::function<void(const ExtTrajectory &)> &maybe,
        OutputFunctional &output) {
        double dist = equalTimeDistance(t, queryTrajectory);
        if (dist < queryDelta) {
            returnTrajectory(t, output);
        } else {
            maybe(t);
        }
    }

    // Query step. The final step for each query is to do a full decision
    // frechet
    // computation.
    // This step contains no additional smart optimization, and so is very slow.
    template <typename OutputFunctional>
    void pruneWithDecisionFrechet(CDFQShortcuts &cdfqs,
                                  const ExtTrajectory &queryTrajectory,
                                  double queryDelta, const ExtTrajectory &t,
                                  squared_distance &dist2,
                                  OutputFunctional &output) {
        if (cdfqs.calculate(queryTrajectory, t, queryDelta, dist2)) {
            returnTrajectory(t, output);
        }
    }

    // -----------------------------------------

    // Only computes reachable part of freespace diagram, uses shortcuts to skip
    // columns in reachable part
    class CDFQShortcuts {
        struct QEntry {
            int start_row_index;
            int end_row_index;
            double lowest_right;
        };

        std::vector<QEntry> queue[2];
        int queueSize[2];

        double computeSegmentFrechet(const Portal &p, int q,
                                     const Vertices &p_array,
                                     const Vertices &q_array,
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
        bool calculate(const Vertices &P, const Vertices &Q, size_t offset_p,
                       size_t offset_q, size_t size_p, size_t size_q,
                       double queryDelta, double baseQueryDelta,
                       const std::map<int, std::vector<Portal>> &portals,
                       squared_distance &dist2) {
            double startDist = dist2(P[offset_p], Q[offset_q]);
            double endDist = dist2(P[size_p - 1], Q[size_q - 1]);
            if (startDist > queryDelta * queryDelta ||
                endDist > queryDelta * queryDelta)
                return false;
            if (size_p <= offset_p + 1 || size_q <= offset_q + 1)
                return false;  // TODO: do we need this? // WM: added offset

            int first = 0;
            int second = 1;

            Range Rf;
            Range Tf;
            Portal choice;

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
            // bool LFree = computeInterval(Q[0], P[0], P[1], queryDelta, Rf);
            queue[first][0].start_row_index = 0;
            queue[first][0].end_row_index = 0;
            queue[first][0].lowest_right = 0;

            queueSize[first] = 1;
            queueSize[second] = 0;

            // For each column
            for (int column = offset_q; column < static_cast<int>(size_q) - 1;
                 column++) {
                if (queueSize[first] == 0) {
                    // nothing reachable anymore
                    return false;
                }
                queueSize[second] = 0;
                int row = queue[first][0].start_row_index;
                int qIndex = 0;
                // while there's reachable cells left in the queue
                while (qIndex < queueSize[first]) {
                    double left_most_top = 2;
                    // start at reachable cell at the head of the queue, and
                    // continue
                    // until
                    // reachability cannot propagate, consuming the queue as we
                    // progress
                    do {
                        // tracks whether we overshoot the queue
                        bool outsideQueue = qIndex >= queueSize[first];
                        // Right edge stored in Rf, RFree = false means not free
                        bool RFree = computeInterval(
                            Q[column + 1], P[row], P[row + 1], queryDelta, Rf);
                        if (RFree) {
                            if (left_most_top <= 1) {
                                double newLR = Rf.start;
                                if (Rf.isComplete() && queueSize[second] > 0 &&
                                    queue[second][queueSize[second] - 1]
                                            .end_row_index == row - 1) {
                                    // complete reachable right means increase
                                    // previous queue
                                    // entry to span this cell
                                    queue[second][queueSize[second] - 1]
                                        .end_row_index = row;
                                } else {
                                    // push to queue
                                    queue[second][queueSize[second]]
                                        .start_row_index = row;
                                    queue[second][queueSize[second]]
                                        .end_row_index = row;
                                    queue[second][queueSize[second]]
                                        .lowest_right = newLR;
                                    queueSize[second]++;
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
                                        if (Rf.isComplete() && newLR == 0.0 &&
                                            queueSize[second] > 0 &&
                                            queue[second][queueSize[second] - 1]
                                                    .end_row_index == row - 1) {
                                            // complete reachable right means
                                            // increase previous queue
                                            // entry to span this cell
                                            queue[second][queueSize[second] - 1]
                                                .end_row_index = row;
                                        } else {
                                            // push to queue
                                            queue[second][queueSize[second]]
                                                .start_row_index = row;
                                            queue[second][queueSize[second]]
                                                .end_row_index = row;
                                            queue[second][queueSize[second]]
                                                .lowest_right = newLR;
                                            queueSize[second]++;
                                        }
                                    }
                                }
                            }
                        }
                        // Top edge stored in Tf, TFree = false means not free
                        bool TFree =
                            computeInterval(P[row + 1], Q[column],
                                            Q[column + 1], queryDelta, Tf);
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
                        if (!outsideQueue && queueSize[second] > 0 &&
                            queue[second][queueSize[second] - 1]
                                    .end_row_index == row &&
                            Rf.end == 1.0) {
                            // jump-off point possible
                            // check if minimum jump distance is big enough
                            int gapSize = queue[first][qIndex].end_row_index -
                                          queue[first][qIndex].start_row_index;
                            if (gapSize > 1) {
                                const auto iter = portals.find(row);
                                if (iter != portals.cend()) {
                                    const std::vector<Portal> &ports =
                                        iter->second;
                                    choice.source = -1;
                                    for (const Portal &p : ports) {
                                        // int jumpSize = p.destination -
                                        // p.source;
                                        // check if jump within range
                                        if (p.destination <=
                                            queue[first][qIndex]
                                                .end_row_index) {
                                            // check if jump distance fits
                                            double segmentFrechet =
                                                computeSegmentFrechet(
                                                    p, column, P, Q, dist2);
                                            if (segmentFrechet + p.distance <=
                                                baseQueryDelta) {
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
                                        queue[second][queueSize[second] - 1]
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

            int endIndex = queueSize[first] - 1;
            if (endIndex == -1) return false;
            bool exit = queue[first][endIndex].start_row_index ==
                            static_cast<int>(size_p) - 2 &&
                        queue[first][endIndex].lowest_right <= 1;
            return exit || (queue[first][endIndex].end_row_index ==
                                static_cast<int>(size_p) - 2 &&
                            queue[first][endIndex].start_row_index !=
                                static_cast<int>(size_p) - 2);
        }

        bool calculate(const ExtTrajectory &P, const ExtTrajectory &Q,
                       double queryDelta, double baseQueryDelta,
                       squared_distance &dist2) {
            return calculate(P.vertices, Q.vertices, 0, 0, P.size(), Q.size(),
                             queryDelta, baseQueryDelta, P.simpPortals, dist2);
        }

        bool calculate(const ExtTrajectory &P, const ExtTrajectory &Q,
                       double queryDelta, squared_distance &dist2) {
            return calculate(P.vertices, Q.vertices, 0, 0, P.size(), Q.size(),
                             queryDelta, queryDelta, P.simpPortals, dist2);
        }
    };

    // Does binary search on integer range (lowerbound, upperbound),
    // accepts lambda function returning whether the given search index
    // satisfies the search criterion.

    void binaryDoubleSearch(const std::function<int(double)> &f,
                            double upperbound, double lowerbound) {
        double rangeLength = upperbound - lowerbound;
        double avg = lowerbound + (rangeLength) / 2;
        int result = f(avg);
        if (result == 1) {
            binaryDoubleSearch(f, upperbound, avg);
        } else if (result == 0) {
            binaryDoubleSearch(f, avg, lowerbound);
        } else {
            return;
        }
    }

    // Does double & search on integer range (lowerbound, upperbound),
    // accepts lambda function returning whether the given search index
    // satisfies the search criterion.
    static int doubleNsearch(const std::function<bool(int)> &f, int start,
                             int end, int doubleNSearchBase,
                             double doubleNSearchExponentStep) {
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
                k = binaryIntSearch(f, k, prevk);
                return k;
            } else {
                if (k == end - 1) {
                    return k;
                }
                prevk = k;
                k += (int)floor(pow(doubleNSearchBase,
                                    doubleNSearchExponentStep * iteration));
                iteration++;
            }
        }
    }

    static int binaryIntSearch(const std::function<bool(int)> &f,
                               int upperbound, int lowerbound) {
        int rangeLength = upperbound - lowerbound;
        if (rangeLength <= 1) {
            return lowerbound;
        }
        int middle = lowerbound + (rangeLength) / 2;
        bool result = f(middle);
        if (result) {
            return binaryIntSearch(f, upperbound, middle);
        } else {
            return binaryIntSearch(f, middle, lowerbound);
        }
    }

    // This class contains all logic need to compute a simplification
    // of any input trajectory using agarwal simplification with
    // double & search. The algorithm uses EqualTimeDistance.h to
    // satisfy the agarwal constraints.
    class AgarwalSimplification {
       public:
        TrajectorySimplification simplify(ExtTrajectory &t,
                                          double simplificationEpsilon,
                                          squared_distance &dist2) {
            TrajectorySimplification simplified(simplificationEpsilon);

            const Vertices &P = t.vertices;
            simplified.vertices.push_back(P[0]);
            simplified.distances.push_back(0);
            simplified.totals.push_back(0);

            int simpSize = 1;

            int rangeStart = 1;
            int prevk = 0;
            while (true) {
                int k = findLastFrechetMatch(
                    P, simplified.vertices, t.totals, simplified.totals,
                    t.distances, simplified.distances,
                    static_cast<int>(simplified.vertices.size()), rangeStart,
                    static_cast<int>(t.size()), prevk, simplificationEpsilon,
                    dist2);
                simpSize++;
                simplified.vertices[simpSize - 1] = P[k];
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
        int findLastFrechetMatch(const Vertices &P, Vertices &simp,
                                 std::vector<double> &ptotals,
                                 std::vector<double> &simptotals,
                                 std::vector<double> &pdists,
                                 std::vector<double> &simpdists, int simpSize,
                                 int start, int end, int prevk, double epsilon,
                                 squared_distance &dist2) {
            simp.push_back(P[0]);
            simpdists.push_back(0);
            simptotals.push_back(0);
            // Use lambda's to easily double & search the function from (start)
            // to
            // (end)
            constexpr int doubleNSearchBase = 2;
            constexpr double doubleNSearchExponentStep = 1.0;
            return doubleNsearch(
                [&](int index) -> bool {
                    simp[simpSize] = P[index];
                    simpdists[simpSize] =
                        std::sqrt(dist2(simp[simpSize], simp[simpSize - 1]));
                    simptotals[simpSize] =
                        simptotals[simpSize - 1] + simpdists[simpSize];
                    double dist = equalTimeDistance(
                        P, simp, ptotals, simptotals, pdists, simpdists,
                        index + 1, simpSize + 1, prevk, simpSize - 1);
                    return dist <= epsilon;
                },
                start, end, doubleNSearchBase, doubleNSearchExponentStep);
        }
    };

    // This class contains all logic need to compute a simplification
    // of any input trajectory using agarwal simplification with
    // double & search. The algorithm uses EqualTimeDistance.h to
    // satisfy the agarwal constraints.
    class ProgressiveAgarwal {
       public:
        TrajectorySimplification simplify(ExtTrajectory &parent,
                                          ExtTrajectory &sourceTrajectory,
                                          double simplificationEpsilon,
                                          squared_distance &dist2) {
            TrajectorySimplification simplified(simplificationEpsilon);

            const Vertices &P = parent.vertices;
            simplified.vertices.push_back(P[0]);
            simplified.distances.push_back(0);
            simplified.totals.push_back(0);
            simplified.sourceIndex.push_back(0);

            int simpSize = 1;

            int rangeStart = 1;
            int prevk = 0;
            while (true) {
                int k = findLastFrechetMatch(
                    P, simplified.vertices, parent.totals, simplified.totals,
                    parent.distances, simplified.distances,
                    static_cast<int>(simplified.vertices.size()), rangeStart,
                    static_cast<int>(parent.size()), prevk,
                    simplificationEpsilon, dist2, parent.sourceIndex,
                    sourceTrajectory, simplified.portals);
                simpSize++;
                simplified.vertices[simpSize - 1] = P[k];
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
        int findLastFrechetMatch(
            const Vertices &P, Vertices &simp, std::vector<double> &ptotals,
            std::vector<double> &simptotals, std::vector<double> &pdists,
            std::vector<double> &simpdists, int simpSize, int start, int end,
            int prevk, double epsilon, squared_distance &dist2,
            std::vector<size_t> &parentSourceIndices,
            ExtTrajectory &sourceTrajectory, std::vector<Portal> &portals) {
            simp.push_back(P[0]);
            simpdists.push_back(0);
            simptotals.push_back(0);
            // Use lambda's to easily double & search the function from (start)
            // to
            // (end)
            constexpr int doubleNSearchBase = 2;
            constexpr double doubleNSearchExponentStep = 1.0;
            return doubleNsearch(
                [&](int index) -> bool {
                    simp[simpSize] = P[index];
                    simpdists[simpSize] =
                        std::sqrt(dist2(simp[simpSize], simp[simpSize - 1]));
                    simptotals[simpSize] =
                        simptotals[simpSize - 1] + simpdists[simpSize];
                    size_t end = 0;
                    size_t start = parentSourceIndices[prevk];
                    if (static_cast<size_t>(index + 1) >=
                        parentSourceIndices.size()) {
                        end = sourceTrajectory.size();
                    } else {
                        end = parentSourceIndices[index + 1];
                    }
                    double dist = equalTimeDistance(
                        sourceTrajectory.vertices, simp,
                        sourceTrajectory.totals, simptotals,
                        sourceTrajectory.distances, simpdists,
                        static_cast<int>(end), simpSize + 1,
                        static_cast<int>(start), simpSize - 1);
                    Portal p;
                    p.source = prevk;
                    p.destination = index;
                    p.distance = dist;
                    portals.push_back(p);
                    return dist <= epsilon;
                },
                start, end, doubleNSearchBase, doubleNSearchExponentStep);
        }
    };
};

}  // namespace bddm
}  // namespace detail
}  // namespace frechetrange

#endif  // TUE_INC
