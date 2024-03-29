# robin-surface-mesh
A program to generate a mesh of the ROBIN body - a generic helicopter body shape

![screenshot](doc/img_good.png?raw=true "ROBIN Body")

### Publication
A full write-up was published as a Technical Note: [Robust and Corrected Coefficients for the
ROtor-Body-INteraction Body, AIAA Journal, Vol. 59, No. 10, October 2021](https://arc.aiaa.org/doi/10.2514/1.J060303).

Cibin Joseph created [genROBIN](https://github.com/cibinjoseph/genROBIN), a Python implementation that exports to more file formats than this version.

### Usage
Here is how you can build and run the program on Linux (or OSX with `g++` installed with Homebrew):

    git clone https://github.com/Applied-Scientific-Research/robin-surface-mesh.git
    cd robin-surface-mesh
    make
    ./makerobin 24 16 12 12

The final command can be used to make the geometry that appears in the image above. The order of command-line arguments is:

    ./makerobin [nx fuselage] [nt fuselage] [nx pylon] [nt pylon]

Where `nx` is the number of lengthwise and `nt` is the number of circumferential elements. Output is to two triangle mesh files in OBJ format containing the fuselage and the pylon.

### Notes
The formulas and coefficients available in the literature are incorrect! They will not only generate NaN values but also not create the correct shape. It seems that no authors mention this, despite echoing the same formulae and coefficients as in the original paper. This program contains the correct coefficients and formulae, and it is hoped that this brings some clarity to the field.
Here is [the original report](https://ntrs.nasa.gov/search.jsp?R=19790017844), note equation 5 and the tables on page 8. Some corrections are made in a [subsequent report](https://ntrs.nasa.gov/search.jsp?R=19870008231), but still the formulae and many coefficients are erroneous.

