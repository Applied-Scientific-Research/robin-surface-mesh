/*
 * makerobin.cpp - generate the ROBIN helicopter body model
 *
 * g++ -o makerobin makerobin.cpp -lm
 */

#include <cstdlib>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
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
  else if (x >= 0.8 and x <= 1.018) idx = 5;
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
  if (!denom) { denom = 1;}
  //std::cout << "Numer=" << numer << " Denom=" << denom << " R=" << (numer / std::pow(denom, 1.0/N)) << std::endl;
  return numer / std::pow(denom, 1.0/N); 
}

void create_mesh(const size_t nx, const size_t nt, std::vector<SupEll> hcoeff, std::vector<SupEll> wcoeff,
                 std::vector<SupEll> zcoeff, std::vector<SupEll> ncoeff, const std::string fileName,
                 int (*getSection)(double), const double xBegin, const double xEnd) {
  // Open file to write to
  std::ofstream file;
  file.open(fileName);
  file << "# Vertices\n";

  std::cout << std::endl << "Generating Nodes" << std::endl << std::endl;
  for (size_t ix=0; ix<nx+1; ix++) {
    
    const double xol = ((xEnd - xBegin) * ix / (double)nx) + xBegin;
    const int fusSec = getSection(xol);
    if (fusSec == -1) {
      std::cout << "ERROR: fusSec == " << fusSec << " x == " << xol << std::endl;
      exit(0);
    }
    //std::cout << "x index=" << ix << " with xol=" << xol << " uses station=" << fusSec << std::endl;

    // compute H, W, Z0, and N from xol and the constants
    //std::cout << "H:" << std::endl;
    const double H  = getsuperval(xol, hcoeff[fusSec]);
    //std::cout << "W:" << std::endl;
    const double W  = getsuperval(xol, wcoeff[fusSec]);
    //std::cout << "Z:" << std::endl;
    const double Z0 = getsuperval(xol, zcoeff[fusSec]);
    //std::cout << "N:" << std::endl;
    const double N  = getsuperval(xol, ncoeff[fusSec]);
    //std::cout << "at xol=" << xol << " have " << H << " " << W << " " << Z0 << " " << N << std::endl;

    for (size_t it=0; it<nt; it++) {
      const double theta = 2.0*3.14159265358979*it/(double)nt;
      // compute r from H, W, N, theta
      const double r = getRadialCoord(H, W, theta, N);
      // compute yol, zol from r, theta, Z0
      //std::cout << "r: " << r << std::endl;
      const double yol = r * std::sin(theta);
      const double zol = r * std::cos(theta) + Z0;
      std::cout /*<< r << "  " */<< xol << " " << yol << " " << zol << std::endl;
      
      // Write vertice to file
      file << "v " << xol << " " << yol << " " << zol << "\n";
      if ((xol == 0) || (xol == 2)) { break; }
    }
  }
  file.close();
}

// execution starts here
int main(int argc, char const *argv[]) {
  std::cout << "Generating ROBIN model" << std::endl;

  // first 4 rows of each section are the fuselage
  //   next 2 more rows each for the pylon

  // fixes:
  // 1) if there's a 0.0 in the second col, then change the 4th and 5th cols to 1.0
  // 2) if there's a 0.0 in C7, change C8 to 1.0, same as above, to prevent nan/inf
  // 3) the 0.4..0.8 section (row 2) coefficients in C1 needed to go into C6
  // 4) C4 is wrong in the first section of fuse and pyl - it needed to be negative

  //std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, 0.0, 0.25, 1.8},
  //                               {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
  //                               {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
  //                               {1.0, -1.0, -0.8, 0.4, 3.0, 0.0, 0.145, 3.0},
  //                               {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.145, 2.0} };
  std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, -0.4, 1.8, 0.0, 0.25, 1.8},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
                                 {1.0, -1.0, -0.8, -0.4, 3.0, 0.0, 0.145, 3.0},
                                 {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.145, 2.0} };

  //std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, 0.4, 2.0, 0.0, 0.25, 2.0},
  //                               {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
  //                               {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
  //                               {1.0, -1.0, -0.8, -0.4, 3.0, 0.0, 0.166, 3.0},
  //                               {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.166, 2.0} };
  std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, -0.4, 2.0, 0.0, 0.25, 2.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
                                 {1.0, -1.0, -0.8, 0.4, 3.0, 0.0, 0.166, 3.0},
                                 {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.166, 2.0} };

  //std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, -0.08, 0.08, 1.8},
  //                               {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.04, -0.04, 0.6},
  //                               {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {1.0, -1.0, -0.8, 1.1, 1.5, 0.065, 0.06, 0.6} };
  std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, -0.4, 1.8, -0.08, 0.08, 1.8},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.04, -0.04, 0.6},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 0.04, 0.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 0.125, 0.0, 1.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.065, 0.06, 0.6} };

  //std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
  //                               {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 0.0, 0.0},
  //                               {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
  //                               {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };
  std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0},
                                 {5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 1.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0},
                                 {0.0, 0.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0} };

  //double tx[3] = {0.5, 0.5, 0.5};
  const size_t nx = 10;
  const size_t nt = 10;
  const std::string fileName = "robin_fuselage.obj";
  const double fusBegin = 0.0;
  const double fusEnd = 2.0;
  const double pylBegin = 0.4;
  const double pylEnd = 1.018;
  const std::string pFileName = "robin_pylon.obj";

  create_mesh(nx, nt, hcoeff, wcoeff, zcoeff, ncoeff, fileName, get_fuselage_section, fusBegin, fusEnd);
  // generate a second closed tri mesh for the pylon, then CSG them together

  std::cout << std::endl << "Generating Triangles" << std::endl;
  // Label faces
  std::ofstream file;
  file.open(fileName, std::ios_base::app);
  file << "\n# Faces\n";
  for (size_t i=2; i<nt+1; i++) {
    file << "f " << 1 << " " << i << " " << i+1 << "\n";
  }
  file << "f " << 1 << " " << nt+1 << " " << 2 << "\n";

  // Link the two rings together to make faces
  for (size_t r=0; r<nx-2; r++) {
    file << "f " << 3+r*nt << " " << 2+r*nt << " " << 2+(r+1)*nt << "\n";
    // Runs through all nodes in the ring
    for (size_t n=1; n<nt; n++) {
      file << "f " << (2+n)+r*nt << " " << (1+n)+(r+1)*nt << " " << (2+n)+(r+1)*nt << "\n";
      file << "f " << (3+n)+r*nt << " " << (2+n)+r*nt << " " << (2+n)+(r+1)*nt << "\n";
    }
    file << "f " << 1+(r+1)*nt << " " << 1+(r+2)*nt << " " << (r+2)*nt << "\n";
    //file << "f " << 2+r*nt << " " << 1+(r+2)*nt << " " << 1+(r+1)*nt << "\n";
    file << "f " << 2+(r+1)*nt << " " << 2+r*nt << " " << 1+(r+1)*nt << "\n";
  }

  size_t lastVert = nt*(nx-1)+2;
  for (size_t i=2; i<nt+1; i++) {
    file << "f " << (i+1)+nt*(nx-2) << " " << i+nt*(nx-2) << " " << lastVert << "\n";
  }
  file << "f " << 2+nt*(nx-2) << " " << lastVert-1 << " " << lastVert << "\n";
  // Done writing faces
  file.close();

  create_mesh(nx, nt, hcoeff, wcoeff, zcoeff, ncoeff, pFileName, get_pylon_section, pylBegin, pylEnd);
  return 0;
}
