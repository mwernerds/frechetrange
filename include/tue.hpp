#ifndef TUE_INC
#define TUE_INC


#include <vector>
#include <thread>
#include <chrono>
#include <mutex>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cfloat>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include <unordered_set>



namespace tue_details{

/*

    This is the collected and reordered (but not redactionally completed) implementation of 

1. Tom van Diggelen    t.w.t.v.diggelen@student.tue.nl
2. Yago Diez           contact@yagodiez.com
3. Kevin Buchin        k.a.buchin@tue.nl
4. Wouter Meulemans    w.meulemans@tue.nl

submitted to the ACM SIGSPATIAL GIS Cup 2017

Refactoring done by Martin in order to have a single-header implementation

*/
// Program-wide defines
#define WRITE_OUTPUT_TO_QUERY true	// true -> outputfiles (result-XXXXX.txt) are created when queries are solved
#define USE_GPU false				// true -> OpenCL is used (DONT FLIP THIS, DOESNT WORK (yet))
#define USE_MULTITHREAD true		// false -> use only one thread, useful to debug concurrency issues
#define USE_FAST_IO true			// true -> file loading is faster, but less robust
#define ONLY_TOTAL_TIMES false		// true -> print diagnostic information
#define USE_FOPEN_S false			// true -> using windows file API

int numSimplifications = 4;
int avgs[4];
int count = 0;
int slotsPerDimension = 500;
double tolerance = 0.00001;

double avgsBBRatio[4];


/* DATA STRUCTURES */
struct Vertex {
	double x;
	double y;

	int trajectoryNumber;
	bool isStart;
};

// Represents a boundingbox around a trajectory or group of trajectories
class BoundingBox {
public:
	double minx = DBL_MAX;
	double miny = DBL_MAX;

	double maxx = -DBL_MAX;
	double maxy = -DBL_MAX;

	void addPoint(double x, double y) {
		if (x < minx) minx = x;
		if (y < miny) miny = y;

		if (x > maxx) maxx = x;
		if (y > maxy) maxy = y;
	}

	double getDiagonal() {
		double width = maxx - minx;
		double height = maxy - miny;

		return sqrt(width*width + height*height);
	}
};


struct Portal {
	int source;
	int destination;
	double distance;
};
class TrajectorySimplification;

class Trajectory {
public: 
	std::string name;
	
	std::vector<Vertex> vertices;
	std::vector<double> distances;
	std::vector<double> totals;
	std::vector<int> sourceIndex;
	std::map<int, std::vector<Portal>> simpPortals;

	int size;
	int uniqueIDInDataset;
	double totalLength;

	BoundingBox *boundingBox;
	std::vector<TrajectorySimplification*> simplifications;

	~Trajectory() {
		delete boundingBox;
	}

	void print() {

		std::cout << "Trajectory: " << name << "\n";
		for (int i = 0; i < size; i++) {
			Vertex v = vertices[i];
			std::cout << i << " " << v.x << " " << v.y << " " << distances[i] << " " << totals[i] << "\n";
		}
	}
};

class TrajectorySimplification : public Trajectory {

public:
	Trajectory* source;
	std::vector<Portal> portals;
	double simplificationEpsilon;
};

bool portalCompare(Portal &lhs, Portal &rhs) { return lhs.destination < rhs.destination; };







// Implements Equal Time Distance algorithm between two trajectories.
// The ETD algorithm computes an approximation of frechet distance by
// taking the 'dog leash' length when traversing two trajectories at 
// the same speed. Used by agarwal and simplification step.
// If found relevant, this function can be optimized using SIMD instructions
static double equalTimeDistance(
	std::vector<Vertex> &pverts, std::vector<Vertex> &qverts,
	std::vector<double> &ptotals, std::vector<double> &qtotals,
	std::vector<double> &pdistances, std::vector<double> &qdistances,
	int psize, int qsize,
	int pstart, int qstart) {

	double pdistOffset = ptotals[pstart];
	double qdistOffset = qtotals[qstart];
	double pdist = ptotals[psize - 1] - pdistOffset;
	double qdist = qtotals[qsize - 1] - qdistOffset;
	double pscale = qdist / pdist;
	int p_ptr = pstart + 1;
	int q_ptr = qstart + 1;

	//startpoints
	double dx = pverts[pstart].x - qverts[qstart].x;
	double dy = pverts[pstart].y - qverts[qstart].y;
	double smax = dx*dx + dy*dy;
	dx = pverts[psize - 1].x - qverts[qsize - 1].x;
	dy = pverts[psize - 1].y - qverts[qsize - 1].y;
	double emax = dx*dx + dy*dy;
	if (qdist == 0 || pdist == 0) return sqrt(std::max(emax, smax));
	double position = 0; // from 0 to 1

	Vertex p_pt;
	Vertex q_pt;

	while (
		!(
			p_ptr == psize - 1
			&&
			q_ptr == qsize - 1
			)
		) {
		double posP = position * pdist;
		double posQ = position * qdist;
		double nextDistP = ptotals[p_ptr] - pdistOffset - posP;
		double nextDistQ = qtotals[q_ptr] - qdistOffset - posQ;

		if (p_ptr == psize - 1) nextDistP = DBL_MAX;
		if (q_ptr == qsize - 1) nextDistQ = DBL_MAX;



		if (nextDistP * pscale < nextDistQ) { // treat P first
			p_pt.x = pverts[p_ptr].x;
			p_pt.y = pverts[p_ptr].y;
			position = (ptotals[p_ptr] - pdistOffset) / pdist;
			double scale = (position * qdist - (qtotals[q_ptr - 1] - qdistOffset)) / qdistances[q_ptr];
			dx = qverts[q_ptr].x - qverts[q_ptr - 1].x;
			dy = qverts[q_ptr].y - qverts[q_ptr - 1].y;
			q_pt.x = qverts[q_ptr - 1].x + dx * scale;
			q_pt.y = qverts[q_ptr - 1].y + dy * scale;
			p_ptr++;
		}
		else { // treat Q first
			q_pt.x = qverts[q_ptr].x;
			q_pt.y = qverts[q_ptr].y;
			position = (qtotals[q_ptr] - qdistOffset) / qdist;
			double scale = (position * pdist - (ptotals[p_ptr - 1] - pdistOffset)) / pdistances[p_ptr];
			dx = pverts[p_ptr].x - pverts[p_ptr - 1].x;
			dy = pverts[p_ptr].y - pverts[p_ptr - 1].y;
			p_pt.x = pverts[p_ptr - 1].x + dx * scale;
			p_pt.y = pverts[p_ptr - 1].y + dy * scale;
			q_ptr++;
		}
		dx = p_pt.x - q_pt.x;
		dy = p_pt.y - q_pt.y;
		double nm = dx*dx + dy*dy;
		if (nm > smax) {
			smax = nm;
		}
	}

	//endpoints
	dx = pverts[p_ptr].x - qverts[q_ptr].x;
	dy = pverts[p_ptr].y - qverts[q_ptr].y;
	double nm = dx*dx + dy*dy;
	if (nm > smax) {
		smax = nm;
	}


	return sqrt(smax);
}

double equalTimeDistance(Trajectory &p, Trajectory &q) {
	return equalTimeDistance(p.vertices, q.vertices, p.totals, q.totals, p.distances, q.distances, p.size, q.size, 0, 0);
}



struct Range {
	double start;
	double end;
};

Range emptyRange = { 0,0 };
Range dontCare = { 0, 0 };

bool isEmpty(Range &r) {
	return r.start == r.end;
}

bool isComplete(Range &r) {
	return r.start == 0.0 && r.end == 1.0;
}

void setRange(Range &r, Range &s) {
	r.start = s.start;
	r.end = s.end;
}

void setRange(Range &r, double start, double end) {
	r.start = start;
	r.end = end;
}

double distSQ(Vertex &p, Vertex &q) {
	double dx = p.x - q.x;
	double dy = p.y - q.y;
	return dx*dx + dy*dy;
}

double clamp01(double input) {
	return input > 1 ? 1 : (input < 0 ? 0 : input);
}

inline bool computeInterval(Vertex &a, Vertex &b1, Vertex &b2, double eps, Range &r) {
	// compute the interval along b1-b2 that is within eps-range of a

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

	// pull out some of the identical computations in here
	double b2m1x = b2.x - b1.x;
	double b2m1y = b2.y - b1.y;
	double b1max = b1.x - a.x;
	double b1may = b1.y - a.y;

	double A = b2m1x * b2m1x + b2m1y * b2m1y;
	double B = 2 * ((b2m1x) * (b1max) + (b2m1y) * (b1may));
	double C = b1max*b1max + b1may * b1may - eps*eps;

	double D = B * B - 4 * A * C;
	if (D < 0) {
		// no solution
		return false;
	}
	else {
		// pull out the Sqrt(D)
		double sqrtD = sqrt(D);
		double t1 = (-B + sqrtD) / (2 * A);
		double t2 = (-B - sqrtD) / (2 * A);
		// rather than doing a swap, if may be faster to check the sign of A before doing the assignment, OR just use min/max, etc 
		double tempt1 = t1;
		t1 = std::min(t1, t2);
		t2 = std::max(tempt1, t2);

		if (t2 < 0 || t1 > 1) {
			return false;
		}
		else {
			r.start = std::max(0.0, t1);
			r.end = std::min(1.0, t2);
			return true;
		}
	}
};

// Adapted from implementation of Yago Diez, could be sped up but is not a bottleneck
class DiHash
{
public:

	int slotsPerDimension; //number of equally spaced subdivisions in each dimension
	double limits[2][2]; // two rows (x,y) and two columns(min,max), keeps the information on limits in each dimension, for 3D, just add one row.

	double tol; // tolerance to prevent numerical representation errors

	std::vector<Vertex>** elements;

	DiHash(BoundingBox &boundingBox, int numC, double iTol) {

		slotsPerDimension = numC;
		tol = iTol;

		// First, locate numerichal limits

		// initialize limits at a three-component two doubles (max and min) vector

		// set initial values for extreme values
		limits[0][0] = boundingBox.minx; //min in X
		limits[0][1] = boundingBox.maxx; //max in X
		limits[1][0] = boundingBox.miny; //min in Y
		limits[1][1] = boundingBox.maxy; //max in Y

		// Now that we have extreme dimensions, set points in its corresponding part of the grid

		elements = new std::vector<Vertex>*[slotsPerDimension];
		for (int i = 0; i < slotsPerDimension; i++) {
			elements[i] = new std::vector<Vertex>[slotsPerDimension];
		}
	}


	void addPoint(Vertex pActual) {
		int x, y;

		x = findSlot(pActual.x, 'x', false);
		y = findSlot(pActual.y, 'y', false);

		elements[x][y].push_back(pActual);

	}

	int findSlot(double val, char type, bool allowOverflow) {
		double min, max;
		int retorn;

		//cout<<"DiHash::findSlot limits "<<endl<<"x: ("<<limits[0][0]<<" , "<<limits[0][1]<<")"<<endl;
		//cout<<"y: ("<<limits[1][0]<<" , "<<limits[1][1]<<")"<<endl;

		switch (type) {
		case 'x':
			//cout<<"DiHash::findSlot x "<<endl;
			min = limits[0][0];
			max = limits[0][1];
			break;
		case 'y':
			//cout<<"DiHash::findSlot y"<<endl;
			min = limits[1][0];
			max = limits[1][1];
			break;
		default:
			std::cout << "DiHash::findSlot(double val, char type) wrong slot ";
			break;
		}
		// check that we are not exactly at the maximum or minimum value
		if (fabs(max - val) < tol) retorn = slotsPerDimension - 1;
		else if (fabs(min - val) < tol) retorn = 0;
		else {
			double pas = (fabs(max - min) / slotsPerDimension);

			//cout<<"DiHash::findSlot computing "normal" return value Math.fabs(max-min)"<<Math.fabs(max-min)<<" pas "<<pas<<endl;

			retorn = (int)((val - min) / pas);
		}

		//cout<<"DiHash::findSlot return valure before checking strange cases (if asked to): "<<retorn<<endl;

		if ((retorn >= slotsPerDimension) || (retorn < 0)) {
			if (!allowOverflow) {
				std::cout << "DiHash::findSlot(double val, char type) wrong return value ";
			}
			else    // in this case, used for searches, we set values out of the values outside the extreme to the last slot
			{
				if (retorn >= slotsPerDimension) retorn = slotsPerDimension - 1;
				else retorn = 0;
			}
		}


		return retorn;

	}
	// Given a type of search (x or y) and a slot number c, return slot c's lower bound
	double slotLowerBound(int c, char type) {
		double min, max;

		switch (type) {
		case 'x':
			min = limits[0][0];
			max = limits[0][1];
			break;
		case 'y':
			min = limits[1][0];
			max = limits[1][1];
			break;
		default:
			std::cout << "DiHash::slotLowerBound wrong slot";
		}

		double pas = (fabs(max - min) / slotsPerDimension);

		return (min + pas * c);
	}

	// Given a type of search (x or y) and a search range (min, max(, return the two indexes of the first and last slots affected by the search
	int* slotsTouched(double min, double max, char type) {
		// first position of the following vector is for the minimum slot affected and second for the maximum
		int* retorn = new int[2];

		retorn[0] = findSlot(min, type, true);
		retorn[1] = findSlot(max, type, true);

		return retorn;
	}

	// Return neighbors at range strictly less than distance for a given point, does not return the query point if present
	void neighbors(Vertex &p, double eps, std::unordered_set<int> &retorn) {
		int* limitsX = slotsTouched(p.x - eps, p.x + eps, 'x');
		int* limitsY = slotsTouched(p.y - eps, p.y + eps, 'y');

		for (int i = limitsX[0]; i <= limitsX[1]; i++) {
			for (int j = limitsY[0]; j <= limitsY[1]; j++) {
				for (Vertex &pActual : elements[i][j]) {
					if (pActual.isStart == p.isStart) {
						double dx = p.x - pActual.x;
						double dy = p.y - pActual.y;
						double distSQ = dx * dx + dy * dy;
						if (distSQ < eps * eps) {
							retorn.insert(pActual.trajectoryNumber);
						}
					}
				}

			}
		}
		delete[] limitsX;
		delete[] limitsY;
	}

	// Return neighbors at range strictly less than distance for a given point, does not return the query point if present
	void neighborsWithCallback(Vertex &p, Vertex& end, double eps, std::vector<Trajectory*> &trajectories, const std::function< void(Trajectory*) >& emit) {
		int* limitsX = slotsTouched(p.x - eps, p.x + eps, 'x');
		int* limitsY = slotsTouched(p.y - eps, p.y + eps, 'y');

		for (int i = limitsX[0]; i <= limitsX[1]; i++) {
			for (int j = limitsY[0]; j <= limitsY[1]; j++) {
				for (Vertex &pActual : elements[i][j]) {
					if (pActual.isStart == p.isStart) {
						double dx = p.x - pActual.x;
						double dy = p.y - pActual.y;
						double distSQ = dx * dx + dy * dy;
						if (distSQ < eps * eps) {
							Trajectory *t = trajectories[pActual.trajectoryNumber];
							Vertex &pend = t->vertices[t->size - 1];
							double dex = end.x - pend.x;
							double dey = end.y - pend.y;
							double disteSQ = dex * dex + dey * dey;
							if (disteSQ < eps * eps) {
								emit(t);
							}
						}
					}
				}

			}
		}
		delete[] limitsX;
		delete[] limitsY;
	}
	


	~DiHash();
};






// Represents a query in a queryset file
// The queryNumber is the index in the query file,
// and should be used when outputting results
struct Query {
	std::string queryTrajectoryFilename;
	double queryDelta;
	int queryNumber;
};


// Only computes reachable part of freespace diagram
class CDFQueued {

private:
	struct QEntry {
		int row_index;
		double lowest_right;
	};

	std::vector<QEntry> queue[2];
	int queueSize[2];

	double dist(Vertex p, Vertex q) {
		double dx = p.x - q.x;
		double dy = p.y - q.y;
		return sqrt(dx*dx + dy*dy);
	}


public:

	int numRows = 0;
	bool calculate(
		std::vector<Vertex> &P, std::vector<Vertex> &Q,
		int offset_p, int offset_q,
		int size_p, int size_q,
		double queryDelta
	) {

		double startDist = dist(P[offset_p], Q[offset_q]);
		double endDist = dist(P[size_p - 1], Q[size_q - 1]);
		if (startDist > queryDelta || endDist > queryDelta) return false;
		if (size_p <= offset_p + 1 || size_q <= offset_q + 1) return false;//TODO: do we need this? // WM: added offset

		int first = 0;
		int second = 1;

		Range Rf;
		Range Tf;

		// ensure queue capacity
		// WM: added offsets
		int max = size_p - offset_p;
		if (size_q - offset_q > max) max = size_q - offset_q;
		if (queue[0].size() < max) {
			for (int i = queue[0].size(); i < max; i++) {
				queue[0].push_back({ 0,0 });
				queue[1].push_back({ 0,0 });
			}
		}

		// setup
		//WM: this is always free space! by check on startDist
		//bool LFree = computeInterval(Q[0], P[0], P[1], queryDelta, Rf); 
		queue[first][0].row_index = 0;
		queue[first][0].lowest_right = 0;

		queueSize[first] = 1;
		queueSize[second] = 0;


		// For each column
		for (int column = offset_q; column < size_q - 1; column++) {
			if (queueSize[first] == 0) {
				// nothing reachable anymore
				return false;
			}
			queueSize[second] = 0;
			int row = queue[first][0].row_index;
			int qIndex = 0;
			// while there's reachable cells left in the queue
			while (qIndex < queueSize[first]) {
				double left_most_top = 2;
				// start at reachable cell at the head of the queue, and continue until
				// reachability cannot propagate, consuming the queue as we progress
				do {
					// tracks whether we overshoot the queue
					bool outsideQueue = qIndex >= queueSize[first];
					// Right edge stored in Rf, RFree = false means not free
					bool RFree = computeInterval(Q[column + 1], P[row], P[row + 1], queryDelta, Rf);
					if (RFree) {
						if (left_most_top <= 1) {
							// push to queue
							queue[second][queueSize[second]].row_index = row;
							queue[second][queueSize[second]].lowest_right = Rf.start;
							queueSize[second]++;
						}
						else {
							// WM: think you should be checking row here as well
							if (!outsideQueue && row == queue[first][qIndex].row_index && queue[first][qIndex].lowest_right <= Rf.end) {
								// push to queue
								queue[second][queueSize[second]].row_index = row;
								queue[second][queueSize[second]].lowest_right = std::max(queue[first][qIndex].lowest_right, Rf.start);
								queueSize[second]++;
							}
						}
					}
					// Top edge stored in Tf, TFree = false means not free
					bool TFree = computeInterval(P[row + 1], Q[column], Q[column + 1], queryDelta, Tf);
					if (!outsideQueue && row == queue[first][qIndex].row_index) {
						// consume the first queue
						qIndex++;
						if (TFree) {
							left_most_top = Tf.start;
						}
						else {
							left_most_top = 2;
						}
					}
					else if (TFree && left_most_top <= Tf.end) {
						left_most_top = std::max(left_most_top, Tf.start);
					}
					else {
						left_most_top = 2;
					}
					// propagated reachability by one cell, so look at next row
					row++;
					numRows++;
				} while (left_most_top <= 1 && row < size_p - 1);
			}

			// swap first and second column
			int temp = first;
			first = second;
			second = temp;
		}

		int endIndex = queueSize[first] - 1;
		return queue[first][endIndex].row_index == size_p - 2 && queue[first][endIndex].lowest_right <= 1;
	}

	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta);
	}
};



// Only computes reachable part of freespace diagram, uses shortcuts to skip columns in reachable part
class CDFQShortcuts {

private:
	struct QEntry {
		int start_row_index;
		int end_row_index;
		double lowest_right;
	};

	std::vector<QEntry> queue[2];
	int queueSize[2];

	double dist(Vertex p, Vertex q) {
		double dx = p.x - q.x;
		double dy = p.y - q.y;
		return sqrt(dx*dx + dy*dy);
	}

	inline double computeSegmentFrechet(
		Portal &p,
		int q,
		std::vector<Vertex> &p_array, std::vector<Vertex> &q_array
	) {
		Vertex &pstart = p_array[p.source];
		Vertex &pend = p_array[p.destination];

		Vertex &qstart = q_array[q];
		Vertex &qend = q_array[q];

		return computeSegmentFrechet(pstart, pend, qstart, qend);
	}

	inline double computeSegmentFrechet(
		Vertex &pstart, Vertex &pend,
		Vertex &qstart, Vertex &qend) {
		double startdx = pstart.x - qstart.x;
		double startdy = pstart.y - qstart.y;
		double enddx = pend.x - qend.x;
		double enddy = pend.y - qend.y;
		double startdist = startdx*startdx + startdy*startdy;
		double enddist = enddx*enddx + enddy*enddy;
		return sqrt(std::max(startdist, enddist));
	}



public:

	int numRows = 0;
	bool calculate(
		std::vector<Vertex> &P, std::vector<Vertex> &Q,
		int offset_p, int offset_q,
		int size_p, int size_q,
		double queryDelta,
		double baseQueryDelta,
		std::map<int, std::vector<Portal>> &portals
	) {
		double startDist = dist(P[offset_p], Q[offset_q]);
		double endDist = dist(P[size_p - 1], Q[size_q - 1]);
		if (startDist > queryDelta || endDist > queryDelta) return false;
		if (size_p <= offset_p + 1 || size_q <= offset_q + 1) return false;//TODO: do we need this? // WM: added offset

		int first = 0;
		int second = 1;

		Range Rf;
		Range Tf;
		Portal choice;

		// ensure queue capacity
		// WM: added offsets
		int max = size_p - offset_p;
		if (size_q - offset_q > max) max = size_q - offset_q;
		if (queue[0].size() < max) {
			for (int i = queue[0].size(); i < max; i++) {
				queue[0].push_back({ 0,0 });
				queue[1].push_back({ 0,0 });
			}
		}

		// setup
		//WM: this is always free space! by check on startDist
		//bool LFree = computeInterval(Q[0], P[0], P[1], queryDelta, Rf); 
		queue[first][0].start_row_index = 0;
		queue[first][0].end_row_index = 0;
		queue[first][0].lowest_right = 0;


		queueSize[first] = 1;
		queueSize[second] = 0;


		// For each column
		for (int column = offset_q; column < size_q - 1; column++) {
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
				// start at reachable cell at the head of the queue, and continue until
				// reachability cannot propagate, consuming the queue as we progress
				do {
					// tracks whether we overshoot the queue
					bool outsideQueue = qIndex >= queueSize[first];
					// Right edge stored in Rf, RFree = false means not free
					bool RFree = computeInterval(Q[column + 1], P[row], P[row + 1], queryDelta, Rf);
					if (RFree) {
						if (left_most_top <= 1) {
							double newLR = Rf.start;
							if (isComplete(Rf) && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row - 1) {
								// complete reachable right means increase previous queue entry to span this cell
								queue[second][queueSize[second] - 1].end_row_index = row;
							}
							else {
								// push to queue
								queue[second][queueSize[second]].start_row_index = row;
								queue[second][queueSize[second]].end_row_index = row;
								queue[second][queueSize[second]].lowest_right = newLR;
								queueSize[second]++;
							}
						}
						else {
							// WM: think you should be checking row here as well
							if (!outsideQueue && row >= queue[first][qIndex].start_row_index && row <= queue[first][qIndex].end_row_index) {
								if (!(row == queue[first][qIndex].start_row_index && queue[first][qIndex].lowest_right > Rf.end)) {
									double prevR = row == queue[first][qIndex].start_row_index ? queue[first][qIndex].lowest_right : 0.0;
									double newLR = std::max(prevR, Rf.start);
									if (isComplete(Rf) && newLR == 0.0 && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row - 1) {
										// complete reachable right means increase previous queue entry to span this cell
										queue[second][queueSize[second] - 1].end_row_index = row;
									}
									else {
										// push to queue
										queue[second][queueSize[second]].start_row_index = row;
										queue[second][queueSize[second]].end_row_index = row;
										queue[second][queueSize[second]].lowest_right = newLR;
										queueSize[second]++;
									}
								}
							}
						}
					}
					// Top edge stored in Tf, TFree = false means not free
					bool TFree = computeInterval(P[row + 1], Q[column], Q[column + 1], queryDelta, Tf);
					if (!outsideQueue && row <= queue[first][qIndex].end_row_index && row >= queue[first][qIndex].start_row_index) {
						if (row == queue[first][qIndex].end_row_index) {
							// consume the first queue
							qIndex++;
						}
						if (TFree) {
							left_most_top = Tf.start;
						}
						else {
							left_most_top = 2;
						}
					}
					else if (TFree && left_most_top <= Tf.end) {
						left_most_top = std::max(left_most_top, Tf.start);
					}
					else {
						left_most_top = 2;
					}
					//try and jump
					if (!outsideQueue && queueSize[second] > 0 && queue[second][queueSize[second] - 1].end_row_index == row && Rf.end == 1) {
						// jump-off point possible
						// check if minimum jump distance is big enough
						int gapSize = queue[first][qIndex].end_row_index - queue[first][qIndex].start_row_index;
						if (gapSize > 1) {
							std::vector<Portal> &ports = portals[row];
							choice.source = -1;
							for (Portal &p : ports) {
								int jumpSize = p.destination - p.source;
								// check if jump within range
								if (p.destination <= queue[first][qIndex].end_row_index) {
									// check if jump distance fits
									double segmentFrechet = computeSegmentFrechet(p, column, P, Q);
									if (segmentFrechet + p.distance <= baseQueryDelta) {
										choice = p;
									}
								}
								else {
									break;
								}
							}
							// JUMP!
							if (choice.source != -1) {
								row = choice.destination - 1;// - 1 to counter ++ later
								queue[second][queueSize[second] - 1].end_row_index = row;
							}
						}
					}
					// propagated reachability by one cell, so look at next row
					row++;
					numRows++;
				} while (left_most_top <= 1 && row < size_p - 1);
			}

			// swap first and second column
			int temp = first;
			first = second;
			second = temp;
		}

		int endIndex = queueSize[first] - 1;
		if (endIndex == -1) return false;
		bool exit = queue[first][endIndex].start_row_index == size_p - 2 && queue[first][endIndex].lowest_right <= 1;
		return exit || (queue[first][endIndex].end_row_index == size_p - 2 && queue[first][endIndex].start_row_index != size_p - 2);
	}

	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta, double baseQueryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta, baseQueryDelta, P.simpPortals);
	}

	bool calculate(Trajectory &P, Trajectory &Q, double queryDelta) {
		return calculate(P.vertices, Q.vertices, 0, 0, P.size, Q.size, queryDelta, queryDelta, P.simpPortals);
	}
};



class FileIO {
	std::vector<Vertex> vertexBuffer;
	std::vector<double> distanceBuffer;
	std::vector<double> totalBuffer;
	std::vector<int> sourceIndex;


public:

	// Parses trajectory file, also computes trajectory metrics
	// TODO: move temporary buffer allocation out of function
	Trajectory* parseTrajectoryFileStreams(std::string filename, int trajectoryNumber) {
		std::ifstream infile(filename);



		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}
		std::string line;
		std::getline(infile, line);
		double x, y, z, w;
		vertexBuffer.clear();
		distanceBuffer.clear();
		totalBuffer.clear();
		sourceIndex.clear();

		distanceBuffer.push_back(0);
		totalBuffer.push_back(0);

		BoundingBox *b = new BoundingBox();

		bool start = true;

		Trajectory *t = new Trajectory();
		t->name = filename;
		Vertex v;

		while (infile >> x >> y >> z >> w)
		{
			v.x = x;
			v.y = y;
			v.trajectoryNumber = trajectoryNumber;
			v.isStart = start;
			start = false;
			//update boundingbox
			b->addPoint(v.x, v.y);
			if (vertexBuffer.empty()) {
				vertexBuffer.push_back(v);
				sourceIndex.push_back(sourceIndex.size());
			}
			else {
				// ignore duplicate verts, they are annoying
				int prevIndex = vertexBuffer.size() - 1;
				Vertex &prev = vertexBuffer[prevIndex];
				if (prev.x != v.x || prev.y != v.y) {
					double dx = v.x - prev.x;
					double dy = v.y - prev.y;
					double dist = sqrt(dx*dx + dy*dy);
					distanceBuffer.push_back(dist);
					totalBuffer.push_back(totalBuffer[prevIndex] + dist);
					vertexBuffer.push_back(v);
					sourceIndex.push_back(sourceIndex.size());
				}
			}
		}
		infile.close();

		t->size = vertexBuffer.size();
		t->totalLength = totalBuffer[t->size - 1];
		t->boundingBox = b;

		t->vertices = vertexBuffer;
		t->distances = distanceBuffer;
		t->totals = totalBuffer;
		t->sourceIndex = sourceIndex;



		return t;
	}

	#define BUFFER_SIZE (1024 * 1024 * 1024)//1Mb
	#define READ_SIZE (1024 * 1024 * 5)//5kb

	//std::vector<char> buffer = std::vector<char>(BUFFER_SIZE);
	char* buffer = nullptr;
	Trajectory* parseTrajectoryFileFast(std::string filename, int trajectoryNumber) {
		if (buffer == nullptr) {
			buffer = new char[BUFFER_SIZE];
		}
#if USE_FOPEN_S
		FILE* file;
		int err = fopen_s(&file, filename.c_str(), "rb");
		if (err != 0) {
			std::cout << "Failed to open (fopens): " << filename << "\n";
			exit(1);
		}
#else
		FILE* file = fopen(filename.c_str(), "rb");
		if (file == NULL) {
			std::cout << "Failed to open (fopen): " << filename << "\n";
			exit(1);
		}
#endif


		double x, y, z, w;
		vertexBuffer.clear();
		distanceBuffer.clear();
		totalBuffer.clear();
		sourceIndex.clear();

		distanceBuffer.push_back(0);
		totalBuffer.push_back(0);

		BoundingBox *b = new BoundingBox();

		bool start = true;

		Trajectory *t = new Trajectory();
		t->name = filename;
		Vertex v;

		int line = 0;

		int offset = 0;
		int reads = 0;

		while (true)
		{
			size_t readLength = fread(&buffer[READ_SIZE*reads], 1, READ_SIZE, file);
			reads++;
			while (true) {
				if (offset == readLength) break;
				line++;
				int eol = -1;
				for (int i = offset; i < offset + readLength; i++) {
					if (buffer[i] == '\n') {
						eol = i;
						break;
					}
				}
				int prevOffset = offset;
				int nextOffset = eol + 1;
				if (eol == -1) break;
				if (line == 1) {
					offset = eol + 1;
					continue;
				}
				// vertex from 0 to nl
				int space = -1;
				for (int i = offset; i < eol; i++) {
					if (buffer[i] == ' ') {
						space = i;
						break;
					}
				}
				offset = nextOffset;
				if (space == -1 || space == prevOffset) break;// bad line
				char* pEnd;
				double x = strtod(&buffer[prevOffset], &pEnd);
				double y = strtod(pEnd, NULL);
				v.x = x;
				v.y = y;
				v.trajectoryNumber = trajectoryNumber;
				v.isStart = start; start = false;
				//update boundingbox
				b->addPoint(v.x, v.y);
				if (vertexBuffer.empty()) {
					vertexBuffer.push_back(v);
					sourceIndex.push_back(sourceIndex.size());
				}
				else {
					// ignore duplicate verts, they are annoying
					int prevIndex = vertexBuffer.size() - 1;
					Vertex &prev = vertexBuffer[prevIndex];
					if (prev.x != v.x || prev.y != v.y) {
						double dx = v.x - prev.x;
						double dy = v.y - prev.y;
						double dist = sqrt(dx*dx + dy*dy);
						distanceBuffer.push_back(dist);
						totalBuffer.push_back(totalBuffer[prevIndex] + dist);
						vertexBuffer.push_back(v);
						sourceIndex.push_back(sourceIndex.size());
					}
				}
			}

			//EOF
			if (readLength != BUFFER_SIZE) {
				break;
			}
		}
		fclose(file);


		t->size = vertexBuffer.size();
		t->totalLength = totalBuffer[t->size - 1];
		t->boundingBox = b;

		t->vertices = vertexBuffer;
		t->distances = distanceBuffer;
		t->totals = totalBuffer;
		t->sourceIndex = sourceIndex;



		return t;
	}

	// delegating function for file loading
	Trajectory* parseTrajectoryFile(std::string filename, int trajectoryNumber) {
		#if USE_FAST_IO
			//return parseTrajectoryFileFast("C:\\frechetdata\\tdrive\\" + filename, trajectoryNumber);
			return parseTrajectoryFileFast(filename, trajectoryNumber);
		#else
			return parseTrajectoryFileStreams(filename, trajectoryNumber);
		#endif
	}


	// Parses query file, does not load query trajectories
	std::vector<Query>* parseQueryFile(char* filename) {
		std::ifstream infile(filename);

		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		double eps;
		std::string queryTrajectoryFileName;

		std::vector<Query> *queries = new std::vector<Query>;

		int queryNumber = 0;
		while (infile >> queryTrajectoryFileName >> eps) {
			Query q;
			q.queryNumber = queryNumber;
			q.queryDelta = eps;
			q.queryTrajectoryFilename = queryTrajectoryFileName;

			queries->push_back(q);
			queryNumber++;
		}
		infile.close();

		return queries;

	}

	// Parses dataset file, does not translate filenames based on dataset location
	// This behavior is consistent with the conversations on the mailing list
	std::vector<std::string>* parseDatasetFile(std::string filename) {
		std::ifstream infile(filename);

		if (!infile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		std::string trajectoryFileName;
		std::vector<std::string> * trajectoryNames = new std::vector<std::string>();

		int trajectoryNumber = 0;
		while (infile >> trajectoryFileName) {
			trajectoryNames->push_back(trajectoryFileName);
		}
		infile.close();

		return trajectoryNames;
	}

	// Writes result trajectories for one query to a file called result-XXXXX.txt 
	void writeQueryOutputFile(Query &q, std::vector<std::string> results) {
		std::ostringstream stringStream;
		stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
		std::string filename = stringStream.str();
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}

		for (std::string trajectoryName : results) {
			outfile << trajectoryName << "\n";
		}
		outfile.close();
	}

	void writeQueryOutputFile(Query &q, std::string results) {
		std::ostringstream stringStream;
		stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
		std::string filename = stringStream.str();
		std::ofstream outfile(filename);

		if (!outfile.is_open()) {
			std::cout << "Failed to open: " << filename << "\n";
			exit(1);
		}
		outfile << results;
		outfile.close();
	}
};// Contains various geometric utility functions used by frechet algorithms.






// Does binary search on integer range (lowerbound, upperbound),
// accepts lambda function returning whether the given search index
// satisfies the search criterion.
int binaryIntSearch(const std::function< bool(int) >& f, int upperbound, int lowerbound) {
	int rangeLength = upperbound - lowerbound;
	if (rangeLength <= 1) {
		return lowerbound;
	}
	int middle = lowerbound + (rangeLength) / 2;
	bool result = f(middle);
	if (result) {
		return binaryIntSearch(f, upperbound, middle);
	}
	else {
		return binaryIntSearch(f, middle, lowerbound);
	}
}

void binaryDoubleSearch(const std::function< int(double) >& f, double upperbound, double lowerbound) {
	double rangeLength = upperbound - lowerbound;
	double avg = lowerbound + (rangeLength) / 2;
	int result = f(avg);
	if (result == 1) {
		binaryDoubleSearch(f, upperbound, avg);
	}
	else if (result == 0){
		binaryDoubleSearch(f, avg, lowerbound);
	}
	else {
		return;
	}
}

// Does double & search on integer range (lowerbound, upperbound),
// accepts lambda function returning whether the given search index
// satisfies the search criterion.
int doubleNsearch(const std::function< bool(int) >& f, int start, int end, int doubleNSearchBase, double doubleNSearchExponentStep) {
	int k = start;
	int prevk = start;
	int iteration = 0;
	while (true) {
		//double
		if (k > end - 1) {
			k = end - 1;
		}
		bool epsValid = f(k);
		if (!epsValid) {
			//binary search
			k = binaryIntSearch(
				f,
				k,
				prevk
				);
			return k;
		}
		else {
			if (k == end - 1) {
				return k;
			}
			prevk = k;
			k += (int)floor(pow(doubleNSearchBase, doubleNSearchExponentStep*iteration));
			iteration++;
		}
	}
}


// This class contains all logic need to compute a simplification
// of any input trajectory using agarwal simplification with
// double & search. The algorithm uses EqualTimeDistance.h to
// satisfy the agarwal constraints.
class AgarwalSimplification {
	int doubleNSearchBase = 2;
	double doubleNSearchExponentStep = 1;

	std::vector<Vertex> simpBuffer;
	std::vector<double> simpDistances;
	std::vector<double> simpTotals;

public:
	TrajectorySimplification* simplify(Trajectory &t, double simplificationEpsilon) {
		TrajectorySimplification* simplified = new TrajectorySimplification();
		simplified->name = t.name + "[simplified]";
		simplified->simplificationEpsilon = simplificationEpsilon;
		simplified->source = &t;

		simplify(t, *simplified, simplificationEpsilon);

		return simplified;
	}

	void simplify(Trajectory &t, TrajectorySimplification &simplification, double simplificationEpsilon) {
		simpBuffer.clear();
		simpDistances.clear();
		simpTotals.clear();

		std::vector<Vertex> &P = t.vertices;
		simpBuffer.push_back(P[0]);
		simpDistances.push_back(0);
		simpTotals.push_back(0);

		int simpSize = 1;

		int rangeStart = 1;
		int prevk = 0;
		while (true) {
			int k = findLastFrechetMatch(P, simpBuffer, t.totals, simpTotals, t.distances, simpDistances, simpSize, rangeStart, t.size, prevk, simplificationEpsilon, simplification.portals);
			simpSize++;
			simpBuffer[simpSize - 1] = P[k];
			if (k == t.size - 1) {
				break;
			}
			prevk = k;
			rangeStart = k + 1;
		}

		simplification.size = simpBuffer.size();


		simplification.vertices = simpBuffer;
		simplification.distances = simpDistances;
		simplification.totals = simpTotals;
	}

private:
	// Finds index k of last vertex v that still satisfies 
	// the simplification epsilon 
	int findLastFrechetMatch(
		std::vector<Vertex> &P, std::vector<Vertex> &simp,
		std::vector<double> &ptotals, std::vector<double> &simptotals,
		std::vector<double> &pdists, std::vector<double> &simpdists,
		int simpSize,
		int start, int end, int prevk, double epsilon,
		std::vector<Portal> &portals) {
		simp.push_back(P[0]);
		simpdists.push_back(0);
		simpTotals.push_back(0);
		// Use lambda's to easily double & search the function from (start) to (end)
		return doubleNsearch(
			[&](int index) -> bool {
				simp[simpSize] = P[index];
				double dx = simp[simpSize].x - simp[simpSize - 1].x;
				double dy = simp[simpSize].y - simp[simpSize - 1].y;
				double d = sqrt(dx*dx + dy*dy);
				simpdists[simpSize] = d;
				simptotals[simpSize] = simptotals[simpSize - 1] + d;
				double dist = equalTimeDistance(
					P, simp,
					ptotals, simptotals,
					pdists, simpdists,
					index + 1, simpSize + 1,
					prevk, simpSize - 1
				);
				return dist <= epsilon;
			},
			start,
			end,
			doubleNSearchBase,
			doubleNSearchExponentStep
		);
	}


};





// This class contains all logic need to compute a simplification
// of any input trajectory using agarwal simplification with
// double & search. The algorithm uses EqualTimeDistance.h to
// satisfy the agarwal constraints.
class ProgressiveAgarwal {
	int doubleNSearchBase = 2;
	double doubleNSearchExponentStep = 1;

	std::vector<Vertex> simpBuffer;
	std::vector<double> simpDistances;
	std::vector<double> simpTotals;
	std::vector<int> sourceIndex;

public:
	TrajectorySimplification* simplify(Trajectory &parent, Trajectory &sourceTrajectory, double simplificationEpsilon) {
		TrajectorySimplification* simplified = new TrajectorySimplification();
		simplified->name = parent.name + "[simplified]";
		simplified->simplificationEpsilon = simplificationEpsilon;
		simplified->source = &parent;

		simplify(parent, *simplified, sourceTrajectory, simplificationEpsilon);
		
		return simplified;
	}

	void simplify(Trajectory &parent, TrajectorySimplification &simplification, Trajectory &sourceTrajectory, double simplificationEpsilon) {
		simpBuffer.clear();
		simpDistances.clear();
		simpTotals.clear();
		sourceIndex.clear();

		std::vector<Vertex> &P = parent.vertices;
		simpBuffer.push_back(P[0]);
		simpDistances.push_back(0);
		simpTotals.push_back(0);
		sourceIndex.push_back(0);

		int simpSize = 1;

		int rangeStart = 1;
		int prevk = 0;
		while (true) {
			int k = findLastFrechetMatch(P, simpBuffer, parent.totals, simpTotals, parent.distances, simpDistances, simpSize, rangeStart, parent.size, prevk, simplificationEpsilon, parent.sourceIndex, sourceTrajectory, simplification.portals);
			simpSize++;
			simpBuffer[simpSize - 1] = P[k];
			sourceIndex.push_back(parent.sourceIndex[k]);
			if (k == parent.size - 1) {
				break;
			}
			prevk = k;
			rangeStart = k + 1;
		}

		simplification.size = simpBuffer.size();


		simplification.vertices = simpBuffer;
		simplification.distances = simpDistances;
		simplification.totals = simpTotals;
		simplification.sourceIndex = sourceIndex;
	}

private:
	// Finds index k of last vertex v that still satisfies 
	// the simplification epsilon 
	int findLastFrechetMatch(
		std::vector<Vertex> &P, std::vector<Vertex> &simp,
		std::vector<double> &ptotals, std::vector<double> &simptotals,
		std::vector<double> &pdists, std::vector<double> &simpdists,
		int simpSize,
		int start, int end, int prevk, double epsilon,
		std::vector<int> &parentSourceIndices,
		Trajectory &sourceTrajectory,
		std::vector<Portal> &portals) {
		simp.push_back(P[0]);
		simpdists.push_back(0);
		simpTotals.push_back(0);
		// Use lambda's to easily double & search the function from (start) to (end)
		return doubleNsearch(
			[&](int index) -> bool {
				simp[simpSize] = P[index];
				double dx = simp[simpSize].x - simp[simpSize - 1].x;
				double dy = simp[simpSize].y - simp[simpSize - 1].y;
				double d = sqrt(dx*dx + dy*dy);
				simpdists[simpSize] = d;
				simptotals[simpSize] = simptotals[simpSize - 1] + d;
				int end = 0;
				int start = parentSourceIndices[prevk];
				if (index + 1 >= parentSourceIndices.size()) {
					end = sourceTrajectory.size;
				}
				else {
					end = parentSourceIndices[index + 1];
				}
				double dist = equalTimeDistance(
					sourceTrajectory.vertices, simp,
					sourceTrajectory.totals, simptotals,
					sourceTrajectory.distances, simpdists,
					end, simpSize + 1,
					start, simpSize - 1
				);
				Portal p;
				p.source = prevk;
				p.destination = index;
				p.distance = dist;
				portals.push_back(p);
				return dist <= epsilon;
			},
			start,
			end,
			doubleNSearchBase,
			doubleNSearchExponentStep
			);
	}


};// Contains the main structure of the algorithm
// The actual logic used to solve each query is in AlgoSteps.h


// Data needed by each worker thread executing the algorithm
// Each thread has seperate algorithm objects so they can use
// their own buffers to store intermediate results.
struct AlgorithmObjects {
	std::ostringstream results;
	std::vector<Trajectory*> candidates;
	BoundingBox bbox;

	FileIO fio;
	AgarwalSimplification agarwal;
	ProgressiveAgarwal agarwalProg;
	CDFQueued cdfq;
	CDFQShortcuts cdfqs;

};

////////// GLOBALS (REMOVE)
// Calculates numSimplification trajectory simplifications for one trajectory
// TODO: improve the epsilon to apply to more types of trajectories,
// TODO: or track the simplification coarseness and adapt epsilon based on it.
void makeSimplificationsForTrajectory(Trajectory &t, double diagonal, AlgorithmObjects &algo, int size) {
	double targets[4] = {.07, .19, .24, .32};

	int targetCounts[4];
	for (int i = 0; i < 4; i++) {
		targetCounts[i] = t.size * targets[i];
		targetCounts[i] = std::max(20, targetCounts[i]);
	}
	targetCounts[0] = std::min(18, targetCounts[0]);//start simple in case dihash is useless

	double diag = t.boundingBox->getDiagonal();
	double lowerBound = diagonal / 100000;
	double upperBound = diagonal / 2;
	int numIterations = 10;

	for (int i = 0; i < size; i++) {
		double targetVertexPercentage = targets[i];
		TrajectorySimplification* simp = nullptr;
		int tries = 0;
		double newUpperbound = 0;
		binaryDoubleSearch(
			[&](double value) -> int {
				newUpperbound = value;
				if (simp != nullptr) delete simp;
				simp = algo.agarwal.simplify(t, value);
				tries++;
				if (tries == 10) {
					return -1;
				}
				else {
					return simp->size > targetCounts[i];
				}
			},
			upperBound,
			lowerBound
		);
		upperBound = newUpperbound;
		numIterations -= 2;
		double ratio = simp->size/(double)t.size;
		avgsBBRatio[i] += newUpperbound / diagonal;
		t.simplifications.push_back(simp);
	}
	count++;

	/*
	for (int i = 0; i < size; i++) {
		TrajectorySimplification* ts = algo.agarwal.simplify(t, simplificationEpsilon / (6 * (i + .25)));
		int diff = ts->size - t.simplifications[i]->size;
		avgs[i] += diff;
	}
	double avg = avgs[2] / (double)count;
	std::cout << "avg " << avg << "\n";
	*/

	// compile portals
	for (int i = 0; i < size; i++) {
		for (Portal &p : t.simplifications[i]->portals) {
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
	//sort portals from small to large
	for (std::map<int, std::vector<Portal>>::iterator iter = t.simpPortals.begin(); iter != t.simpPortals.end(); ++iter)
	{	
		std::vector<Portal> &k = iter->second;
		std::sort(k.begin(), k.end(), portalCompare);
	}
}

// Calculates numSimplification trajectory simplifications for one trajectory
// TODO: improve the epsilon to apply to more types of trajectories,
// TODO: or track the simplification coarseness and adapt epsilon based on it.
void makeSourceSimplificationsForTrajectory(Trajectory &t, Trajectory &source, double diagonal, AlgorithmObjects &algo, int size) {
	for (int i = 0; i < size; i++) {
		double eps = diagonal * (avgsBBRatio[i]/count);
		t.simplifications.push_back(algo.agarwalProg.simplify(t, source, eps));
	}
	// compile portals
	for (int i = 0; i < size; i++) {
		for (Portal &p : t.simplifications[i]->portals) {
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
	//sort portals from small to large
	for (std::map<int, std::vector<Portal>>::iterator iter = t.simpPortals.begin(); iter != t.simpPortals.end(); ++iter)
	{
		std::vector<Portal> &k = iter->second;
		std::sort(k.begin(), k.end(), portalCompare);
	}
}




// Datastructure containing worker threads
std::vector<std::thread*> threads;
struct AlgoData {
	std::vector<Query> *queries;
	std::vector<Trajectory*> *trajectories;
	std::vector<std::string> *trajectoryNames;
	int numTrajectories;
	DiHash* diHash;
	FileIO fio;
	BoundingBox* boundingBox;
	volatile int startedSolving = 0;
	volatile int startedSimplifying = 0;
	int numWorkers;
};

// pre-processing steps --------------------------------------------------------------
// Query step. Does rangequeries for start/endpoints of dataset. Adds all found trajectories
// to candidates.
void collectDiHashPoints(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, const std::function< void(Trajectory*) >& emit) {
	Vertex start = queryTrajectory.vertices[0];
	Vertex end = queryTrajectory.vertices[queryTrajectory.size - 1];

	a->diHash->neighborsWithCallback(start, end, q.queryDelta, *a->trajectories, [&](Trajectory* t) -> void{
		emit(t);
	});
}



// Preprocessing step. Inserts start and endpoints in a regular grid so they can be used
// for range queries later.
void addPtsToDiHash(AlgoData &a) {
	a.diHash = new DiHash(*a.boundingBox, slotsPerDimension, tolerance);
	for (Trajectory *t : *a.trajectories) {
		if (t != nullptr) {
			a.diHash->addPoint(t->vertices[0]);
			a.diHash->addPoint(t->vertices[t->size - 1]);
		}
	}
}
 

void makeSimplificationsForTrajectory(Trajectory &t, AlgorithmObjects &algo) {
	makeSimplificationsForTrajectory(t, t.boundingBox->getDiagonal(), algo, numSimplifications);
}

void loadAndSimplifyTrajectory(std::string &tname, int tIndex, AlgorithmObjects &algo, AlgoData &a) {
	Trajectory *t = algo.fio.parseTrajectoryFile(tname, tIndex);
	if (t->size == 1) {
		delete t;
		// ugly but necessary
		a.trajectories->at(tIndex) = nullptr;
		return;
	}
	algo.bbox.addPoint(t->boundingBox->minx, t->boundingBox->miny);
	algo.bbox.addPoint(t->boundingBox->maxx, t->boundingBox->maxy);
	makeSimplificationsForTrajectory(*t, algo);
	a.trajectories->at(tIndex) = t;

}

// Mutex guarding access to the queryset from the worker threads
std::mutex simplificationMTX;
// Number of queries allocated to a worker as one 'job'
int simplificationBatchSize = 20;

// Returns a query index for a worker to solve, locking the query set
int getConcurrentTrajectory(AlgoData *a) {
	simplificationMTX.lock();
	if (a->startedSimplifying > a->trajectories->size()) {
		simplificationMTX.unlock();
		return -1;
	}
	int returnTrajectory = a->startedSimplifying;
	a->startedSimplifying += simplificationBatchSize;
	simplificationMTX.unlock();
	return returnTrajectory;
}

void simplificationWorker(AlgoData *a, AlgorithmObjects *algo) {
	int current = getConcurrentTrajectory(a);
	std::vector<std::string> &trajectories = *a->trajectoryNames;
	while (current != -1) {
		int limit = simplificationBatchSize;
		if (current + simplificationBatchSize > trajectories.size()) {
			limit = trajectories.size() - current;
		}
		for (int step = 0; step < limit; step++) {
			std::string &t = trajectories[current + step];
			loadAndSimplifyTrajectory(t, current + step, *algo, *a);
		}
		current = getConcurrentTrajectory(a);
	}
}


std::vector<std::thread*> simplificationThreads;
std::vector<AlgorithmObjects*> algos;








// SENSIBILITY CHECKED MARKER




// Preprocessing step. Calculates simplifications for all trajectories in the dataset
void constructSimplifications(AlgoData &a) {
	a.trajectories = new std::vector<Trajectory*>();
	a.trajectories->resize(a.numTrajectories);
	for (int i = 0; i < a.numWorkers; i++) {
		AlgorithmObjects *algo = new AlgorithmObjects();
		std::thread *t = new std::thread(simplificationWorker, &a, algo);
		simplificationThreads.push_back(t);
		algos.push_back(algo);
	}
	for (int i = 0; i < a.numWorkers; i++) {
		(*simplificationThreads[i]).join();
		AlgorithmObjects *algo = algos[i];
		a.boundingBox->addPoint(algo->bbox.minx, algo->bbox.miny);
		a.boundingBox->addPoint(algo->bbox.maxx, algo->bbox.maxy);
		delete simplificationThreads[i];
		delete algo;
	}
	algos.clear();
	simplificationThreads.clear();
	std::cout << "done";
}




// All data needed by the algorithm to solve a specific query file
// Also contains structures needed for preprocessing
 
// Included here to avoid include problem
// #include "AlgoSteps.h"
// @NOTE MIGHT BE RELEVANT

// Does all needed preprocessing for the given the dataset
void preprocessDataSet(AlgoData *a) {
	constructSimplifications(*a);
	addPtsToDiHash(*a);
}


/*
/*
    PRUNING STRATEGIES
    ==================================
    Lorem
*/

// Query step. For each trajectory T in the dataset and query trajectory Q, this step
// compares successive simplifications of T and Q with continuous decision frechet.
// Each comparison can result in YES, NO, or MAYBE.
// YES   -> remove from candidates, add to results
// NO    -> remove from candidates
// MAYBE -> try next simplification, if none are left, continue to next algorithm step
void pruneWithSimplifications(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t, 
	const std::function< void(Trajectory*) >& maybe, const std::function< void(Trajectory*) >& result) {
	bool broke = false;
	for (int i = 0; i < numSimplifications; i++) {

		double decisionEpsilonLower = q.queryDelta
			- queryTrajectory.simplifications[i]->simplificationEpsilon
			- t->simplifications[i]->simplificationEpsilon;

		double decisionEpsilonUpper = q.queryDelta
			+ queryTrajectory.simplifications[i]->simplificationEpsilon
			+ t->simplifications[i]->simplificationEpsilon;

		double dist = equalTimeDistance(*t->simplifications[i], *queryTrajectory.simplifications[i]);

		if (dist < decisionEpsilonLower) {
			result(t);
			broke = true;
			break;
		}

		if (decisionEpsilonLower > 0) {
			bool r = algo->cdfqs.calculate(*queryTrajectory.simplifications[i], *t->simplifications[i], decisionEpsilonLower, q.queryDelta);
			if (r) {
				result(t);
				broke = true;
				break;
			}
		}
		if (decisionEpsilonUpper > 0) {
			bool r = algo->cdfqs.calculate(*queryTrajectory.simplifications[i], *t->simplifications[i], decisionEpsilonUpper, q.queryDelta);
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

// Query step. Uses equal time distance as an upperbound for the actual frechet distance
// If ETD(P, Q) <= queryDelta then CDF(P,Q) <= queryDelta. With P in dataset and Q query trajectory.
void pruneWithEqualTime(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t,
	const std::function< void(Trajectory*) >& maybe, const std::function< void(Trajectory*) >& result) {
	double dist = equalTimeDistance(*t, queryTrajectory);
	if (dist < q.queryDelta) {
		result(t);
	}
	else {
		maybe(t);
	}
}


// Query step. The final step for each query is to do a full decision frechet computation.
// This step contains no additional smart optimization, and so is very slow.
void pruneWithDecisionFrechet(AlgoData *a, Query &q, AlgorithmObjects *algo, Trajectory &queryTrajectory, Trajectory *t,
	const std::function< void(Trajectory*) >& result) {
	bool r = algo->cdfqs.calculate(queryTrajectory, *t, q.queryDelta);
	if (r) {
		result(t);
	}
}


/*
    End Pruning
    ==================================

*/

// Solves a single query, calls functions in AlgoSteps.h
// Before solving, also loads the query trajectory (since it may 
// not be present in the dataset) and constructs simplifications for it.
void solveQuery(AlgoData *a, Query &q, AlgorithmObjects *algo) {
	Trajectory *queryTrajectory = algo->fio.parseTrajectoryFile(q.queryTrajectoryFilename, -1);


	double diagonal = queryTrajectory->boundingBox->getDiagonal();
	makeSourceSimplificationsForTrajectory(*queryTrajectory, *queryTrajectory, diagonal, *algo, numSimplifications);
	for (int i = 1; i < numSimplifications; i++) {
		makeSourceSimplificationsForTrajectory(*queryTrajectory->simplifications[i], *queryTrajectory, diagonal, *algo, i-1);
	}



#if WRITE_OUTPUT_TO_QUERY 
	std::ostringstream stringStream;
	stringStream << "result-" << std::setfill('0') << std::setw(5) << q.queryNumber << ".txt";
	std::string filename = stringStream.str();
	std::ofstream outfile(filename);

	if (!outfile.is_open()) {
		std::cout << "Failed to open: " << filename << "\n";
		exit(1);
	}
#endif

	int results = 0;
	int dihash = 0;
	int simp = 0;
	int et = 0;

	const std::function< void(Trajectory*) >& result = [&](Trajectory *t) -> void{
		results++;
#if WRITE_OUTPUT_TO_QUERY 
		outfile << t->name << "\n";
#endif
	};

	collectDiHashPoints(a, q, algo, *queryTrajectory, [&](Trajectory *t) -> void {
		dihash++;
		pruneWithSimplifications(a, q, algo, *queryTrajectory, t, [&](Trajectory *t) -> void {
			simp++;
			pruneWithEqualTime(a, q, algo, *queryTrajectory, t, [&](Trajectory *t) -> void {
				et++;
				pruneWithDecisionFrechet(a, q, algo, *queryTrajectory, t, result);
			}, result);
		}, result);
	});

#if WRITE_OUTPUT_TO_QUERY 
	outfile.close();
#endif

	//cleanup queryTrajectory, doesn't work from destructor somehow
	for (int i = 0; i < queryTrajectory->simplifications.size(); i++) {
		Trajectory *s = queryTrajectory->simplifications[i];
		for (int j = 0; j < s->simplifications.size(); j++) {
			delete s->simplifications[j];
		}
		delete queryTrajectory->simplifications[i];
	}
	delete queryTrajectory;

	return;

}

// Mutex guarding access to the queryset from the worker threads
std::mutex queryMtx;
// Number of queries allocated to a worker as one 'job'
int querySteps = 20;

// Returns a query index for a worker to solve, locking the query set
int getConcurrentQuery(AlgoData *a) {
	queryMtx.lock();
	if (a->startedSolving > a->queries->size()) {
		queryMtx.unlock();
		return -1;
	}
	int returnQuery = a->startedSolving;
	a->startedSolving += querySteps;
	queryMtx.unlock();
	if (returnQuery % 100 == 0) {
		std::cout << " --- Solving: " << returnQuery << "\n";
	}
	return returnQuery;
}

// Function executed by the worker threads, obtains a query index
// solves a fixed number of queries from that index, then
// tries to obtain a new index. When no queries are left, it exits
// and merges its statistics with the complete statistics.
void worker(AlgoData *a, AlgorithmObjects *algo) {
	int current = getConcurrentQuery(a);
	std::vector<Query> &queries = *a->queries;
	while (current != -1) {
		int limit = querySteps;
		if (current + querySteps > queries.size()) {
			limit = queries.size() - current;
		}
		for (int step = 0; step < limit; step++) {
			Query &c = queries[current + step];
			solveQuery(a, c, algo);
		}
		current = getConcurrentQuery(a);
	}
	delete algo;
}


void printMS(std::string msg, long ms) {
	double timeSec = ms / 1000.0;
	std::cout << msg << ": " << timeSec << " sec \n";
}

void print(std::string msg, double value) {
	std::cout << msg << ": " << value << "\n";
}



// Spins up all worker threads, waits for them to complete,
// then prints statistics.
void solveQueries(AlgoData *a) {
	for (int i = 0; i < a->numWorkers; i++) {
		std::thread *t = new std::thread(worker, a, new AlgorithmObjects());
		threads.push_back(t);
	}
	for (int i = 0; i < a->numWorkers; i++) {
		(*threads[i]).join();
		delete threads[i];
	}
}

void cleanup(AlgoData *a) {
	//TODO: improve code by moving deallocations here
	//TODO: not strictly necessary because program exits
}

// Entrypoint for the algorithm
void runAlgorithm(AlgoData *a) {
	long timeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	std::cout << " - Preprocess\n";
	preprocessDataSet(a);
	long ptimeMS = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	std::cout << " - Solve\n";
	solveQueries(a);
	std::cout << " - Cleanup\n";
	cleanup(a);
	long total = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1) - timeMS;
	long totalp = std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1) - ptimeMS;
	printMS("TOTAL", total);
	printMS("TOTAL_SOLVE", totalp);

}


// Contains the logic used to solve each query,
// called from the main algorithm at Algorithm.h
// Each step changes the set of candidate and result sets,
// where the candidate set indicates a trajectory that MAY
// become a result, but we aren't sure yet, and the result set 
// contains the actual results.







 









// runtime steps ---------------------------------------------------------------------








#include <algorithm>
#include <vector>
#include <stdio.h>






#include <algorithm>
#include <vector>
#include <stdio.h>






#include <stdio.h>
#include <vector>
#include <unordered_set>




#include <functional>
#include <algorithm>






#include <vector>
#include <math.h>
#include <algorithm>





// Contains wrappers for all file IO needed for the algorithm
// Such as loading trajectories, query files, dataset files, etc



} // namespace tue
#endif //TUE
