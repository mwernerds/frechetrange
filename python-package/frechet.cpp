/*
Frechet Python Package


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


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> point; 
typedef bg::model::box<point> box;
typedef bg::model::polygon<point, false, false> polygon; // ccw, open polygo
typedef bg::model::multi_polygon<polygon> multipolygon;
typedef bg::model::linestring<point> linestring;
typedef std::pair<point, size_t> value;

typedef bgi::rtree< value, bgi::rstar<16, 4> > rtree_t;



namespace bp = boost::python;
namespace np = boost::python::numpy;

// wrap c++ array as numpy array
static boost::python::object wrap(double* data, npy_intp size) {
  using namespace boost::python;

  npy_intp shape[1] = { size }; // array size
  PyObject* obj = PyArray_New(&PyArray_Type, 1, shape, NPY_DOUBLE, // data type
                              NULL, data, // data pointer
                              0, NPY_ARRAY_CARRAY, // NPY_ARRAY_CARRAY_RO for readonly
                              NULL);
  handle<> array( obj );
  return object(array);
}

// wrap 2D c++ array as numpy array
static boost::python::object wrap2D(double* data, npy_intp h, npy_intp w) {
  using namespace boost::python;

  npy_intp shape[2] = { h,w }; // array size
  PyObject* obj = PyArray_New(&PyArray_Type, 2, shape, NPY_DOUBLE, // data type
                              NULL, data, // data pointer
                              0, NPY_ARRAY_CARRAY, // NPY_ARRAY_CARRAY_RO for readonly
                              NULL);
  handle<> array( obj );
  return object(array);
}


template<typename func>
static void map_matrix(np::ndarray &self, func f)
{
   auto nd = self.get_nd();
   if (nd != 2) throw(std::runtime_error("2D array expected"));
   auto s1 = self.strides(0);
   auto s2 = self.strides(1);
   auto data = self.get_data();
//   std::cout << s1 << " stride " << s2 << std::endl;
   for (int i1=0; i1 < self.shape(0); i1++)
   {
     for(int i2=0; i2 < self.shape(1); i2++)
     {
         auto offset = i1 * s1 + i2 * s2;
	 double *d = reinterpret_cast<double *> (data + offset);
	 f(i1,i2,*d);
     }
   }
   

}

//typedef NumericMatrix::ConstRow point_adapter;

typedef point point_adapter;






/*struct trajectory_adapter {
  const NumericMatrix &_m;

  typedef point_adapter value_type;

  trajectory_adapter(const NumericMatrix &m) : _m(m){};
  point_adapter operator[](size_t idx) const { return _m.row(idx); }

  size_t size() const { return _m.nrow(); }
};

//@todo: check for 2D in the R interface

struct get_adapter_coord {
  template <size_t dim> static double get(const point_adapter &c) {
    return c[dim];
  }
};

frechetrange::detail::dv::frechet_distance<2, get_adapter_coord>
    fd;

*/

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

#include "../include/frechetrange/frechetrange.hpp"

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







namespace bp =  boost::python;

BOOST_PYTHON_MODULE(frechet) {
    import_array();
    np::initialize();

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
       });

       

}

