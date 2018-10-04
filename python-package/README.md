# Micro-Manual:

Install python, numpy (pip) and some other stuff.
Download boost from boost.org (>=1.67). boost from many distributions (Debian, Ubuntu, etc.) is too old!


Then, bootstrap like this

    ./bootstrap.sh --exec-prefix=. --with-libraries=python --libdir=/usr/local/lib --includedir=./inc

build like

    ./b2

and install like

    sudo ./b2 install
    ldconfig

Don't forget to adapt the paths in the Makefile. Copy the .so to the directory where your python file is. Nothing else needed.

We are planning to automate the build process completely with pip, but therefore, we need a way to find a compiled version of 
boost::numpy, which is quite new.

