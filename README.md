# robin-surface-mesh
A program to generate a mesh of the ROBIN body - a generic helicopter body shape

### Usage
Here is how you can build and run the program on Linux (or OSX with `g++` installed with Homebrew):

    git clone https://github.com/Applied-Scientific-Research/robin-surface-mesh.git
    cd robin-surface-mesh
    make
    ./makerobin > out.obj

Output is to a triangle mesh file in OBJ format containing both the fuselage and the pylon. You can change the number of lengthwise and circumferential stations on the command-line. See a list of options with:

    ./makerobin -h

### Notes
The formulas and coefficients available in the literature are incorrect! They will not only generate NaN values but also not create the correct shape. It seems that no authors mention this, despite echoing the same formulae and coefficients as in the original paper. This program contains the correct coefficients and formulae, and it is hoped that this brings some clarity to the field.
Here is [the original report](https://ntrs.nasa.gov/search.jsp?R=19790017844), note equation 5 and the tables on page 8.

