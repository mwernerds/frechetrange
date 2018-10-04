MYBOOSTLIB=/usr/local/lib
CFLAGS=-std=c++11 -Ofast -fPIC 

frechet.so: frechet.o
	g++ -shared -o frechet.so frechet.o -lpython2.7 -L$(MYBOOSTLIB) -lboost_python -lboost_numpy27 -lboost_system 
frechet.o: frechet.cpp 
	g++ $(CFLAGS) -o frechet.o -c frechet.cpp  -I/local_home/wern_m3/boost/boost_1_67_0 -I/usr/include/python2.7 -I/usr/lib/python2.7/dist-packages/numpy/core/include/ 
clean:
	rm -f *.o *.so
