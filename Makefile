CXX=g++
CXXFLAGS=-O2 -std=c++17
LIBS=-lm

all : makerobin

makerobin : makerobin.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean : 
	rm -f *.o makerobin
