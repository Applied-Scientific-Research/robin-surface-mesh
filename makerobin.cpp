/*
 * makerobin.cpp - generate the ROBIN helicopter body model
 *
 * g++ -o makerobin makerobin.cpp -lm
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>

int getsection (const double x) {
  int idx = -1;
  if (x >= 0.0 and x < 0.4) idx = 0;
  else if (x >= 0.4 and x < 0.8) idx = 1;
  else if (x >= 0.8 and x < 1.9) idx = 2;
  else if (x >= 1.9 and x <= 2.0) idx = 3;
  else idx = -1;
  return idx;
}

using SupEll = std::array<double,8>;

double getsuperval (const double x, const SupEll& c) {
  //std::cout << "se is " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << std::endl;
  // (x/l + C3)/C4
  double innerB = (x-c[2])/c[3];
  if (innerB < 0.0) {
    innerB = 0;
  } else {
    innerB = pow(innerB, c[4]);
  }
  std::cout << "comps are " << (x+c[2])/c[3] << " " << c[4] << " " << innerB << std::endl;
  // C1 + C2 * [innerB]^ C5
  double outerB = c[0]+c[1]*innerB;
  if (outerB < 0.0) {
    outerB = 0;
  } else {
    outerB = pow(outerB, 1./c[7]);
  }
  return c[5] + c[6]*outerB; // nan
  
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], c[4]), 1./c[7]);	// nan
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], c[4]), 1./2.0);	// nan
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], 2.0), 1./c[7]);	// fine
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x-c[2])/c[3], 1.8), 1./c[7]);		// fine
  //return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], 2.0), 1./2.0);	// fine
}

double getRadialCoord(double H, double W, double theta, double N) {
  double numerator = 0.25*H*W;
  if (numerator < 0.0) {
    numerator = 0.0;
  } else {
    numerator = std::pow(numerator, N);
  }
  double denomH= 0.5*H*std::sin(theta);
  if (denomH < 0.0) {
    denomH = 0;
  } else {
    denomH = std::pow(denomH, N);
  }
  double denomW= 0.5*W*std::cos(theta);
  if (denomH < 0.0) {
    denomW = 0;
  } else {
    denomW = std::pow(denomW, N);
  }

  double denom = denomH + denomW;
  // If the denominator is 0 or both the numerator and denominator are negative
  if ((denom == 0.0) || (numerator * denom < 0.0)) {
   return 0.0;
  } else {
    return std::pow(numerator / denom, 1.0/N);
  }
}

// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << "Generating ROBIN model" << std::endl;

  std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, 0.0, 0.25, 1.8},
                                 {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };

  std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, 0.4, 2.0, 0.0, 0.25, 2.0},
                                 {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };

  std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, -0.08, 0.08, 1.8},
                                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };

  std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
                                 {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };

  //double tx[3] = {0.5, 0.5, 0.5};
  const size_t nx = 100;
  const size_t nt = 80;

  std::cout << std::endl << "generating nodes" << std::endl;
  for (size_t ix=1; ix<=nx ; ++ix) {
    const double xol = 2.0 * ix / (double)nx;
    const int isec = 0;//getsection(xol);
    // compute H, W, Z0, and N from xol and the constants
    const double H  = getsuperval(xol, hcoeff[isec]);
    const double W  = getsuperval(xol, wcoeff[isec]);
    const double Z0 = getsuperval(xol, zcoeff[isec]);
    const double N  = getsuperval(xol, ncoeff[isec]);
    std::cout << "at xol=" << xol << " have " << H << " " << W << " " << Z0 << " " << N << std::endl;
    for (size_t it=0; it<nt ; ++it) {
      const double theta = 2.0*3.14159265358979*it/(double)nt;
      // compute r from H, W, N, theta
      const double r = getRadialCoord(H, W, theta, N);
      /* const double r = std::pow( std::pow(0.25*H*W, N) /
                                (std::pow(0.5*H*std::sin(theta), N)+
                                 std::pow(0.5*W*std::cos(theta), N))
                                , 1.0/N); */
      // compute yol, zol from r, theta, Z0
      std::cout << "r: " << r << std::endl;
      const double yol = r * std::sin(theta);
      const double zol = r * std::cos(theta) + Z0;
      std::cout << "  " << xol << " " << yol << " " << zol << std::endl;
      exit(0);
    }

    exit(0);
    // make the triangles for this band
  }

  std::cout << std::endl << "generating triangles" << std::endl;

  return 0;
}

