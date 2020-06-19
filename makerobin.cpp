/*
 * makerobin.cpp - generate the ROBIN helicopter body model
 *
 * g++ -o makerobin makerobin.cpp -lm
 */

#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <array>

int get_fuselage_section (const double x) {
  int idx = -1;
  if (x >= 0.0 and x < 0.4) idx = 0;
  else if (x >= 0.4 and x < 0.8) idx = 1;
  else if (x >= 0.8 and x < 1.9) idx = 2;
  else if (x >= 1.9 and x <= 2.0) idx = 3;
  return idx;
}

int get_pylon_section (const double x) {
  int idx = -1;
  if (x >= 0.4 and x < 0.8) idx = 4;
  else if (x >= 0.8 and x < 1.018) idx = 5;
  return idx;
}

using SupEll = std::array<double, 8>;

double getsuperval (const double x, const SupEll& c) {
  //std::cout << "se is " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << std::endl;
  //std::cout << "comps are " << ((x+c[2])/c[3]) << " " << c[4] << std::endl;
  double tmp = x+c[2];
  //std::cout << "x+c3=" << tmp;
  tmp = tmp/c[3];
  //std::cout << " /c4=" << tmp;
  tmp = std::pow(tmp, c[4]);
  //std::cout << " ^c5=" << tmp;
  tmp = tmp*c[1];
  //std::cout << " *c2=" << tmp;
  tmp = tmp+c[0];
  //std::cout << " +c1=" << tmp;
  tmp = std::pow(tmp, 1.0/c[7]);
  //std::cout << "^(1/c8)=" << tmp;
  tmp = tmp*c[6];
  //std::cout << " *c7=" << tmp;
  tmp = tmp + c[5];
  //std::cout << " +c6=" << tmp << std::endl;
  return c[5] + c[6]*std::pow(std::max(0.0,c[0]+c[1]*std::pow((x+c[2])/c[3], c[4])), 1.0/c[7]);	// nan
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], c[4]), 1./c[7]);	// nan
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], c[4]), 1./2.0);	// nan
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], 2.0), 1./c[7]);	// fine
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x-c[2])/c[3], 1.8), 1./c[7]);		// fine
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], 2.0), 1./2.0);	// fine
}

double getRadialCoord(double H, double W, double theta, double N) {
  double numer = 0.25*H*W;
  // Note the new std::abs - this is to ensure that values are positive, we really only compute one quadrant
  double denom = std::pow(0.5*H*std::abs(std::sin(theta)), N) + std::pow(0.5*W*std::abs(std::cos(theta)), N);
  //std::cout << "Numer=" << numer << " Denom=" << denom << " R=" << (numer / std::pow(denom, 1.0/N)) << std::endl;
  return numer / std::pow(denom, 1.0/N); 
}
// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << "Generating ROBIN model" << std::endl;

  // first 4 rows are the fuselage
  // need to add 2 more rows each for the pylon

  // fixes:
  // 1) if there's a 0.0 in the second col, then change the 4th and 5th cols to 1.0
  // 2) if there's a 0.0 in C7, change C8 to 1.0, same as above, to prevent nan/inf
  // 3) the 0.4..0.8 section (row 2) coefficients in C1 needed to go into C6
  // 4) C4 is wrong in the first sections - it needed to be negative

  //std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, 0.0, 0.25, 1.8},
  //                               {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
  //                               {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0} };
  std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, -0.4, 1.8, 0.0, 0.25, 1.8},
                                 {0.0, 0.0, 0.0, 1.0, 1.0, 0.25, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0} };

  //std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, 0.4, 2.0, 0.0, 0.25, 2.0},
  //                               {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
  //                               {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0} };
  std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, -0.4, 2.0, 0.0, 0.25, 2.0},
                                 {0.0, 0.0, 0.0, 1.0, 1.0, 0.25, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0} };

  //std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, -0.08, 0.08, 1.8},
  //                               {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.04, -0.04, 0.6},
  //                               {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
  std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, -0.4, 1.8, -0.08, 0.08, 1.8},
                                 {0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.04, -0.04, 0.6},
                                 {0.04, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0} };

  //std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
  //                               {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 0.0, 0.0},
  //                               {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
  std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 1.0, 5.0, 0.0, 1.0},
                                 {5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 1.0, 1.0},
                                 {2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0} };

  //double tx[3] = {0.5, 0.5, 0.5};
  const size_t nx = 40;
  const size_t nt = 40;

  std::cout << std::endl << "generating nodes" << std::endl << std::endl;
  //for (size_t ix=0; ix<nx+1; ix++) {
  for (size_t ix=1; ix<200; ix++) {

    const double xol = 2.0 * ix / (double)nx;
    const int isec = get_fuselage_section(xol);
    if (isec == -1) {
      std::cout << "ERROR: isec == " << isec << std::endl;
      exit(0);
    }
    //std::cout << "x index=" << ix << " with xol=" << xol << " uses station=" << isec << std::endl;

    // compute H, W, Z0, and N from xol and the constants
    //std::cout << "H:" << std::endl;
    const double H  = getsuperval(xol, hcoeff[isec]);
    //std::cout << "W:" << std::endl;
    const double W  = getsuperval(xol, wcoeff[isec]);
    //std::cout << "Z:" << std::endl;
    const double Z0 = getsuperval(xol, zcoeff[isec]);
    //std::cout << "N:" << std::endl;
    const double N  = getsuperval(xol, ncoeff[isec]);
    //std::cout << "at xol=" << xol << " have " << H << " " << W << " " << Z0 << " " << N << std::endl;

    for (size_t it=0; it<nt ; it++) {
      const double theta = 2.0*3.14159265358979*it/(double)nt;
      // compute r from H, W, N, theta
      const double r = getRadialCoord(H, W, theta, N);
      // compute yol, zol from r, theta, Z0
      const double yol = r * std::sin(theta);
      const double zol = r * std::cos(theta) + Z0;
      std::cout /*<< r */<< "  " << xol << " " << yol << " " << zol << std::endl;
      //exit(0);
    }

    //exit(0);
    // make the triangles for this band
  }

  // generate a second closed tri mesh for the pylon, then CSG them together

  std::cout << std::endl << "generating triangles" << std::endl;

  return 0;
}

