CXX=g++
CXXFLAGS=-O2 -std=c++17
LIBS=-lm

all : makerobin FreemanMineck1979 PhelpsBerry1987

% : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean : 
	rm -f *.o makerobin
