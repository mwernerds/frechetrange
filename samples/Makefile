CXX=g++
CXXFLAGS=-std=c++11 -Wall -Werror -pedantic-errors
MULTITHREADED=-DENABLE_MULTITHREADING -lpthread

all:
	$(CXX) -o plain plain.cpp $(CXXFLAGS)
	$(CXX) -o bg bg.cpp $(CXXFLAGS)
	$(CXX) -o grid grid.cpp $(CXXFLAGS) $(MULTITHREADED)
	$(CXX) -o spatial_index spatial_index.cpp $(CXXFLAGS)
	$(CXX) -o tue tue.cpp $(CXXFLAGS) -Wno-sign-compare $(MULTITHREADED)

clean:
	rm -f plain bg grid spatial_index tue

format:
	clang-format -i *.cpp ../include/*.hpp
