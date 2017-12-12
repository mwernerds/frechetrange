
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


typedef frechetrange::detail::duetschvahrenhold::Grid<   
     ContainerAdapter, std::function<double (ContainerAdapterRow, ContainerAdapterRow)>,decltype(getx),decltype(gety)>  dv_grid_t;

std::vector<dv_grid_t > gGrids;
   


// [[Rcpp::export]]
size_t internal_prepareRangeQueryDataset(List &l)
{
    // as far as I understood:    grid(meshSize, squared_dist, getx, gety);
    double meshSize = 1;
    std::function<double (ContainerAdapterRow, ContainerAdapterRow)> d = [](ContainerAdapterRow r1, ContainerAdapterRow r2){return ((r2[0]-r1[0])*(r2[0]-r1[0])+(r2[1]-r1[1])*(r2[1]-r1[1]));};


    dv_grid_t G(meshSize,
     dist_sqr,
     getx,gety);

/*  
    @Fabian
    this does not work, as a lot of deleted constructors in Grid are fired. 

    Can you help me? 
*/

/*
  for (const NumericMatrix &t: l)
    {
	G.insert(ContainerAdapter(t));
    }
  */   
   
    gGrids.push_back(G);
    return (gGrids.size() -1);
}


