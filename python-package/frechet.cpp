/*
Frechet Python Package

TODO: unnötige Abhängigkeiten 
- Dritter Decider ist vermischt mit Index, daher nicht sinnvoll isolierbar
- template parameter DIM, coordinate_getter, squared_distance
- template param niederlande: braucht zusätzlich punkt_typ
*/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include<iostream>

#include<boost/python.hpp>
#include<boost/python/numpy.hpp>
#include<boost/python/numpy/ndarray.hpp>
#include <numpy/ndarrayobject.h> 
#include <numpy/ndarraytypes.h> 

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iomanip>      // std::setprecision
#include <limits>
#include <boost/geometry/algorithms/distance.hpp> 

#include "../include/frechetrange/frechetrange.hpp"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bp =  boost::python;


typedef bg::model::point<double, 2, bg::cs::cartesian> point;

typedef bg::model::box<point> box;
typedef bg::model::polygon<point, false, false> polygon; // ccw, open polygo
typedef bg::model::multi_polygon<polygon> multipolygon;
typedef bg::model::linestring<point> linestring;
typedef std::pair<point, size_t> value;

typedef bgi::rtree< value, bgi::rstar<16, 4> > rtree_t;

namespace bp = boost::python;
namespace np = boost::python::numpy;

typedef point point_adapter;


struct trajectory_adapter {
   const np::ndarray &_m;
   const int s1,s2; // the strides
   const char *data;
   
   typedef point_adapter value_type; // @todo
 
   trajectory_adapter(const np::ndarray &m): _m(m), s1(m.strides(0)), s2(m.strides(1)), data(m.get_data()){
   };
   
   point_adapter operator[](size_t idx) const {
      assert(idx < _m.shape(0));
      double c1 = *const_cast<double *>(reinterpret_cast<const double *>( data + idx * s1));
      double c2 = *const_cast<double *>(reinterpret_cast<const double *>( data + idx * s1 + s2));
      return point(c1,c2);
   }

  size_t size() const { return _m.shape(0); }
};


struct get_adapter_coord {
  template <size_t dim> static double get(const point_adapter &c) {
    return bg::get<dim> (c);
  }
};


using std::cout;
using std::endl;


frechetrange::detail::dv::frechet_distance<2, get_adapter_coord>   fd;
frechetrange::detail::bb::frechet_distance<2, get_adapter_coord>   fd2;

class FrechetDecider
{
    public:
       FrechetDecider()
       {
	    std::cout << "Frechet Decider instantiated" << std::endl;
       }

       bool decide(np::ndarray t1, np::ndarray t2, double dist, std::string decider)
	{
	    if (decider == "duetschvahrenhold") return decide_dv(t1,t2,dist);
	    if (decider == "baldusbringman") return decide_bb(t1,t2,dist);

	    return true;
	}

	bool decide_dv(np::ndarray t1, np::ndarray t2, double dist)
	{
	    trajectory_adapter c1(t1), c2(t2);
	    return fd.is_bounded_by(c1,c2,dist);
 	}
	
	bool decide_bb(np::ndarray t1, np::ndarray t2, double dist)
	{
	    trajectory_adapter c1(t1), c2(t2);
	    return fd2.is_bounded_by(c1,c2,dist);
 	}

	
}; // FrechetDecider



///////////////////// TYPES FOR THE GRID INDICES AND THE COPY API


constexpr size_t g_DIMENSIONS = 2;
template <size_t dims> using _point_t = std::array<double, dims>;

struct get_point_coord {
  template <size_t dim> static double get(const _point_t<g_DIMENSIONS> &p) {
    return std::get<dim>(p);
  }
};

template <size_t dims> using _trajectory_t = std::vector<_point_t<dims>>;

template <size_t dims>
void _copyMatrixToTrajectory(const np::ndarray &m, _trajectory_t<dims> &t) {
  char *row = m.get_data();
  const int s1 = m.strides(0), s2 = m.strides(1);
  for (size_t i = 0; i < m.shape(0); i++) {
    char *p = row;
    for (size_t d = 0; d < dims; d++)
    {
      t[i][d] = *reinterpret_cast< double *> (p);
      p += s2;
    }
    row += s1;
  }
}
template<size_t dims>
void _trajectoryToMatrix(const _trajectory_t<dims> &t, np::ndarray &m) {
  char *row = m.get_data();
  const int s1 = m.strides(0), s2 = m.strides(1);
  for (size_t i = 0; i < m.shape(0); i++) {
    char *p = row;
    for (size_t d = 0; d < dims; d++)
    {
      *reinterpret_cast< double *> (p) = t[i][d];
      p += s2;
    }
    row += s1;
  }
}




/*
template <size_t dims>
NumericMatrix _trajectory_to_matrix(const _trajectory_t<dims> &t) {
  NumericMatrix m(t.size(), dims);
  for (size_t i = 0; i < t.size(); i++)
    for (size_t d = 0; d < dims; d++)
      m(i, d) = t[i][d];
  return m;
}

template <size_t dims> class append_to_list {
  List *_list;

public:
  append_to_list(List &l) : _list(&l) {}

  void operator()(const _trajectory_t<dims> &t) {
    _list->push_back(_trajectory_to_matrix<dims>(t));
  }
};


*/




template <size_t dims> class grid_data {
public:
  typedef frechetrange::detail::dv::grid<
      dims, _trajectory_t<dims>, get_point_coord>
      grid_type;

  grid_data() : _pGrid(nullptr), _trajectories() {}
  ~grid_data() { delete _pGrid; }

  grid_type &grid() { return *_pGrid; }

  size_t add_trajectory(const np::ndarray &m) {
    _trajectories.emplace_back(m.shape(0));
    _copyMatrixToTrajectory(m, _trajectories.back());
    return (_trajectories.size() - 1);
  }

  std::string as_ssv()
  {
    // This is meant mainly for debug, it gets back the data as a new numpy array.
    std::stringstream ss;
    for (size_t i=0; i< _trajectories.size(); i++)
    {
	auto &t =_trajectories[i];
	for (auto &p: t)
	    ss << i << " " << std::get<0>(p) << " " << std::get<1>(p) << std::endl;
    }
     return ss.str();
  }

  size_t size() { return _trajectories.size(); }

  void clear() { _trajectories.clear(); }

  void build_index(double meshSize) {
    _pGrid = new grid_type(meshSize);
    _pGrid->reserve(size());
    for (_trajectory_t<dims> &t : _trajectories)
      _pGrid->insert(std::move(t));
    clear();
    _pGrid->build_index();
  }


  bp::object query(const np::ndarray &m, double dist)
  {
      if (_pGrid== nullptr) throw(std::runtime_error("You need to run build_index before you can query"));
      _trajectory_t<dims> t(m.shape(0));
      _copyMatrixToTrajectory(m, t);
//      _copyMatrixToTrajectory<dims>(m, t);
       std::vector<_trajectory_t<dims>> resultSet;
       auto inserter = [&](const _trajectory_t<dims> &t) {resultSet.push_back(t);};
       grid().range_query(t, dist, inserter);
       bp::list l;
       for (auto &t:resultSet)
       {
	    npy_intp  shape[2] = { static_cast<npy_intp> (t.size()),2 }; // array size
	    np::ndarray ta = np::zeros(bp::make_tuple(t.size(),2),np::dtype::get_builtin<double>());
	    _trajectoryToMatrix(t,ta);
	    l.append(ta);
       }
       return bp::object(l);
  }

  

private:
  std::vector<_trajectory_t<dims>> _trajectories;
  grid_type *_pGrid;
};


typedef grid_data<2> DVGrid;


BOOST_PYTHON_MODULE(frechet) {
    import_array();
    np::initialize();

    bp::class_<DVGrid>("DVGrid")
       .def("add",+[](DVGrid &self,bp::object ot1) {
	   np::ndarray at1 = np::from_object(ot1,np::dtype::get_builtin<double>());
	   self.add_trajectory(at1);
	})
     .def("as_ssv",+[](DVGrid &self)->std::string {
           return self.as_ssv();
	})
     .def("build_index",+[](DVGrid &self, double mesh_size) {
          self.build_index(mesh_size);
	  })
     .def("query",+[](DVGrid &self, bp::object ot1, double distance)->bp::object {
           np::ndarray at1 = np::from_object(ot1,np::dtype::get_builtin<double>());
          return self.query(at1, distance);
	  })
	;

    bp::class_<FrechetDecider>("FrechetDecider")
       .def("decide",+[](FrechetDecider &self,bp::object ot1,bp::object ot2, double d, std::string decider)->bool {
	   np::ndarray at1 = np::from_object(ot1,np::dtype::get_builtin<double>());
	   np::ndarray at2 = np::from_object(ot2,np::dtype::get_builtin<double>());	
	    return self.decide(at1,at2,d,decider);

       })
       .def("decide_dv",+[](FrechetDecider &self, bp::object ot1,bp::object ot2, double d)->bool
       {
	   np::ndarray at1 = np::from_object(ot1,np::dtype::get_builtin<double>());
	   np::ndarray at2 = np::from_object(ot2,np::dtype::get_builtin<double>());	
	    if (at1.get_nd() != 2) throw(std::runtime_error("2D array expected"));
	    if (at2.get_nd() != 2) throw(std::runtime_error("2D array expected"));
	    return self.decide_dv(at1,at2,d);
       })
       .def("decide_bb",+[](FrechetDecider &self, bp::object ot1,bp::object ot2, double d)->bool
       {
	   np::ndarray at1 = np::from_object(ot1,np::dtype::get_builtin<double>());
	   np::ndarray at2 = np::from_object(ot2,np::dtype::get_builtin<double>());	
	    if (at1.get_nd() != 2) throw(std::runtime_error("2D array expected"));
	    if (at2.get_nd() != 2) throw(std::runtime_error("2D array expected"));
	    return self.decide_bb(at1,at2,d);
       })
       

 
       ;
}

