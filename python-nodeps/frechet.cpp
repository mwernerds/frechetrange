#include <stdio.h>
#include <Python.h>


#include "../include/frechetrange/frechetrange.hpp"
#include<iostream>

using std::cout;
using std::endl;

typedef std::pair<double, double> point;
typedef std::vector<point> trajectory;

struct get_coordinate {
    template <size_t dim>
    static double get(const point &p) {
        return std::get<dim>(p);
    }
};

double frechet_distance_impl(const trajectory &t1, const trajectory &t2, double low, double high)
{
    //  trajectory t1 = {{0, 0}, {0, 1}, {0, 2}};
    //  trajectory t2 = {{1, 1}, {2, 2}, {1, 3}};
/*    trajectory t1 = {{0, 0}, {10, 0}};  // straight line along X-axis
    trajectory t2 = {
        {0, 1}, {5, 10}, {10, 1}};  // triangle above. Tip is 10 higher.
*/   // true Frechet: 10
    frechetrange::detail::dv::frechet_distance<
        2, get_coordinate,
        frechetrange::euclidean_distance_sqr<2, get_coordinate>>
        fd;
    while(!fd.is_bounded_by(t1, t2, high)) // make sure frechet distance is somewhere inside
    {
	std::cout << "Extending interval above " << high << std::endl;
	low = high;
	high *=2 ;
    }
    double m;
    // refine estimate the distance by interval cutting
    {
        std::pair<double, double> interval = {low, high};
        while ((interval.second - interval.first) > 1e-6) {
            m = (interval.first + interval.second) / 2;

            if (fd.is_bounded_by(t1, t2, m)) {
                interval.second = m;
            } else {
                interval.first = m;
            }
        }
	      cout << "Final Interval: [" << interval.first << "," << interval.second
	<< "]" << endl;
    }
    std::cout << "Found some value " << m << std::endl;
    return m; 
}

// Module method definitions
static PyObject* frechet_world(PyObject *self, PyObject *args) {
    printf("Frechet, world!\n");
    Py_RETURN_NONE;
}

static PyObject* frechet(PyObject *self, PyObject *args) {
    const char* name;
    if (!PyArg_ParseTuple(args, "s", &name)) {
        return NULL;
    }

    printf("Frechet, %s!\n", name);
    Py_RETURN_NONE;
}

trajectory parse_python_trajectory(PyObject *obj)
{
    trajectory ret;
   PyObject *seq = PySequence_Fast(obj, "expected a sequence");
    int len = PySequence_Size(seq);
    for (int i = 0; i < len; i++) {
	PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
	// this must be an iterator
	PyObject *iter = PyObject_GetIter(item);
	if (!iter) {
	    throw(std::runtime_error("Expecting list of iteratable"));
	}
	PyObject *next = PyIter_Next(iter);
	if (!next)
	    throw(std::runtime_error("Did not find something"));
	if (!PyFloat_Check(next))
	    throw(std::runtime_error("Not a native float, be sure to send floats, not ints"));
	double x= PyFloat_AsDouble(next);
	
	next = PyIter_Next(iter);
	if (!next)
	    throw(std::runtime_error("Did not find something"));
	if (!PyFloat_Check(next))
	    throw(std::runtime_error("Not a native float, be sure to send floats, not ints"));
	double y= PyFloat_AsDouble(next);
	ret.push_back({x,y});

        /* DON'T DECREF item here */
    }
    Py_DECREF(seq);
    return ret;
}



static PyObject* frechet_distance(PyObject *self, PyObject *args) {
    PyObject *obj, *obj2;

    if (!PyArg_ParseTuple(args, "OO", &obj, &obj2)) {
	throw(std::runtime_error("Args must be two 2D iteratable of floats"));
    }
    trajectory t1 = parse_python_trajectory(obj);
    trajectory t2 = parse_python_trajectory(obj2);
    for (const auto &t: t1)
	std::cout << "T1:" << t.first << ";" << t.second << std::endl;

    for (const auto &t: t2)
	std::cout << "T2:" << t.first << ";" << t.second << std::endl;
    
    double d = frechet_distance_impl(t1,t2,0.0, 100.0); // put a sensible upper bound here to control speed
    std::cout << "Result: " << d << std::endl;
    return Py_BuildValue("f",d);
}



// Method definition object for this extension, these argumens mean:
// ml_name: The name of the method
// ml_meth: Function pointer to the method implementation
// ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
// ml_doc:  Contents of this method's docstring
static PyMethodDef frechet_methods[] = { 
    {   
        "frechet_world", frechet_world, METH_NOARGS,
        "Print 'frechet world' from a method defined in a C extension."
    },  
    {   
        "frechet_distance", frechet_distance, METH_VARARGS,
        "Compute Frechet distance."
    },  
    {   
        "frechet", frechet, METH_VARARGS,
        "Print 'frechet xxx' from a method defined in a C extension."
    },  
    {NULL, NULL, 0, NULL}
};

// Module definition
// The arguments of this structure tell Python what to call your extension,
// what it's methods are and where to look for it's method definitions
static struct PyModuleDef frechet_definition = { 
    PyModuleDef_HEAD_INIT,
    "frechet",
    "A Python module that prints 'frechet world' from C code.",
    -1, 
    frechet_methods
};

// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_frechet(void) {
    Py_Initialize();
    return PyModule_Create(&frechet_definition);
}
