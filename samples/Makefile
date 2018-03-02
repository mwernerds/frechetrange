CXX=g++
CXXFLAGS=-std=c++11 -Wall -Werror -pedantic-errors
MULTITHREADED=-DENABLE_MULTITHREADING -lpthread

all:
	$(CXX) -o plain plain.cpp $(CXXFLAGS)
	$(CXX) -o bg bg.cpp $(CXXFLAGS)
	$(CXX) -o grid grid.cpp $(CXXFLAGS) $(MULTITHREADED)
	$(CXX) -o spatial_index spatial_index.cpp $(CXXFLAGS)
	$(CXX) -o spatial_hash spatial_hash.cpp $(CXXFLAGS) -Wno-sign-compare $(MULTITHREADED)

clean:
	rm -f plain bg grid spatial_index spatial_hash

format:
	clang-format -i *.cpp ../include/frechetrange/*.hpp ../include/frechetrange/detail/bb/*.hpp ../include/frechetrange/detail/bddm/*.hpp ../include/frechetrange/detail/dv/*.hpp
