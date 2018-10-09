from setuptools import setup, Extension
import numpy as np;
import os;


boost_lib = os.getenv('BOOST_LIB', '/usr/local/lib')
boost_inc = os.getenv('BOOST_INC','/local_home/wern_m3/boost/boost_1_67_0')
    

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
     name='frechet',  
     version='0.1',
#     scripts=['f'] ,
     author="GIS Cup 2018, Martin Werner",
     author_email="martin@martinwerner.de",
     description="This package contains a python package providing Frechet Range Queries as implemented for the ACM SIGSPATIAL GIS Cup 2017",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/mwernerds/frechetrange.git", 
     #packages=setuptools.find_packages()
     ext_modules=[
        Extension("frechet",
            sources=[
                     'frechet.cpp',
            ],
            include_dirs=[
                os.path.join(np.get_include(), 'numpy'),
                boost_inc
            ],
            library_dirs = [os.getcwd(),boost_lib],  # path to .a or .so file(s)
            extra_compile_args=["-std=c++11"],
            language='c++11',
        ),
    ]
#     classifiers=[
#         "Programming Language :: Python :: 3",
#         "License :: OSI Approved :: MIT License",
#         "Operating System :: OS Independent",
#     ],
 )
