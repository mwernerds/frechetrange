
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include<functional>
#include "../../../include/frechetrange.hpp"
using namespace Rcpp;

struct ContainerAdapterRow{
  const NumericMatrix &m;
  size_t row;
  ContainerAdapterRow(const NumericMatrix &_m, size_t _row): m(_m),row(_row) {};
  double operator[] (size_t col) const
  {
    return m(row,col);
  }
      
};

struct ContainerAdapter{

  const NumericMatrix &m;
  typedef ContainerAdapterRow value_type;
  ContainerAdapter(NumericMatrix &_m): m(_m) {};
  ContainerAdapter(const NumericMatrix &_m): m(_m) {};
  const ContainerAdapterRow operator[](size_t row) const{
  
     return ContainerAdapterRow (m,row);;
     }


  size_t size() const{
       return m.nrow();
  }
};

//@todo: check for 2D in the R interface

auto  dist_sqr = [](ContainerAdapterRow r1, ContainerAdapterRow r2){return ((r2[0]-r1[0])*(r2[0]-r1[0])+(r2[1]-r1[1])*(r2[1]-r1[1]));};
auto getx =	[](ContainerAdapterRow r1){return r1[0];};
auto gety =	[](ContainerAdapterRow r1){return r1[1];};

frechetrange::detail::duetschvahrenhold::FrechetDistance<decltype(dist_sqr),decltype(getx),decltype(gety)>
       fd(dist_sqr,getx,gety);

frechetrange::detail::bringmanbaldus::FrechetDistance<
	ContainerAdapter, ContainerAdapter::value_type, double,
	 decltype(getx),decltype(gety),decltype(dist_sqr)>
	 fd2( getx,gety,dist_sqr);


// [[Rcpp::export]]
bool internal_frechet_decide_dv(NumericMatrix &t1, NumericMatrix &t2, double eps)
{
    ContainerAdapter c1(t1),c2(t2); 
    return fd.isBoundedBy(c1,c2,eps);
}



// [[Rcpp::export]]
bool internal_frechet_decide_bb(NumericMatrix &t1, NumericMatrix &t2, double eps)
{
    ContainerAdapter c1(t1),c2(t2); 
    return fd.isBoundedBy(c1,c2,eps);
}

/// Grid will cache the data outside of R to have valid container as some of the implementations
// are not yet compatible. This is a huge memory-overhead and must be well-documented.
// otoh, it gives us highly efficient storage and allows to control memory layout, e.g. through
// boost::object_pool, where wanted.

template<size_t _DIM=2>
class GridDataset{
   public:
        size_t DIM = _DIM;
	typedef std::vector<double> point;
	typedef std::function<double(point, point)> distance_functional_type;
	typedef std::vector<point> trajectory;

	distance_functional_type squared_dist = [](point p1, point p2) {
  // provide a squared distance functional for the point type
  return (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
};
	std::function<double(point p)> getx = [](point p) { return p[0]; };
	std::function<double(point p)> gety = [](point p) { return p[1]; };


   

       std::vector<trajectory> DB;
       size_t addTrajectory(NumericMatrix m)
       {
	   // @todo: track the distance thresholds min and max on insertion
          trajectory t(m.nrow());
	   for (size_t i=0; i < m.nrow(); i++)
	   {
	       t[i].resize(DIM);
	       for (size_t d = 0; d < DIM; d++)
	          t[i][d] = m(i,d);
	   }
	   DB.push_back(t);
	   return (DB.size() -1);
       }

       size_t size(){
	return DB.size();
       }

       void clear() {
	  DB.clear();
       }

/*
*/
    typedef 	frechetrange::detail::duetschvahrenhold::Grid<
            trajectory,
            distance_functional_type,
            std::function<double(point p)>,
           std::function<double(point p)>> grid_type;
 grid_type *pGrid;

	grid_type &grid()
	{
	   return *pGrid;
	}
 
	void buildIndex()
	{
            double meshSize = 2;
	    pGrid = new grid_type(meshSize, squared_dist, getx, gety);
	    auto &grid = *pGrid;
	    grid.reserve(size()); // not mandatory, but adviced in case of many inserts
	    for (auto &t: DB)
		 grid.insert(t);
	  grid.optimize(); // not mandatory, but adviced after completing all/most inserts

	}

       
};



std::vector<GridDataset<2>> gData;

#define ASSERT_VALID_DATASET(k) {if (k < 0 || k >= gData.size()) throw(std::runtime_error("Invalid Handle. Create one with createGridDataset"));}


// [[Rcpp::export]]
size_t internal_createGridDataset()
{
    gData.push_back(GridDataset<2>());
    return gData.size() -1;
}


// [[Rcpp::export]]
size_t internal_addTrajectory (size_t gds, NumericMatrix &m)
{
   ASSERT_VALID_DATASET(gds);
   return gData[gds].addTrajectory(m);
}

// [[Rcpp::export]]
size_t internal_clearDataset(size_t gds)
{
   ASSERT_VALID_DATASET(gds);
   gData[gds].clear();
}

// [[Rcpp::export]]
bool internal_createIndex_dv(size_t gds)
{
   ASSERT_VALID_DATASET(gds);
   Rcout << "Creating Index from " << gData[gds].size() << " trajectories";
   gData[gds].buildIndex();
}

template<typename traj>
NumericMatrix _trajectory2matrix(traj &t)
{
    // @todo: assert non-empty t
   NumericMatrix m(t.size(), t[0].size());
   for(size_t i=0; i < t.size(); i++)
     for (size_t j=0; j < t[i].size(); j++)
       m(i,j) = t[i][j];
   return m;
}

// [[Rcpp::export]]
List internal_gridRangeQuery(size_t gds, NumericMatrix m, double eps,
                             bool materialize = true)
{
  ASSERT_VALID_DATASET(gds);
  auto &D = gData[gds];  
  decltype(gData)::value_type::trajectory t(m.nrow());
  for (size_t i = 0; i < m.nrow(); i++) {
    t[i].resize(D.DIM);
    for (size_t d = 0; d < D.DIM; d++)
      t[i][d] = m(i, d);
  }
  auto results = D.grid().rangeQuery(t, eps);
  List ret(results.size());
  for (size_t i=0; i < results.size(); i++)
  {
     ret[i] = _trajectory2matrix(*results[i]);    
  }
  return ret;  
}
