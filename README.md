# robin-surface-mesh
A program to generate a mesh of the ROBIN body - a generic helicopter body shape

### usage
This code is incomplete! The formulas and coefficients from the paper generate NaN values. Here is how you can build and run the program:

    g++ -o makerobin makerobin.cpp -lm && ./makerobin

### notes
See [the original paper](http://markjstock.org/transfer/20000057579.pdf), page 11, for the formula. The formula seems wrong.
