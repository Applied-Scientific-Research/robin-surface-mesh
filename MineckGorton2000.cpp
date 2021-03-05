/*
 * MineckGorton2000.cpp - generate the ROBIN helicopter body model using 2000 coefficients and formulae
 *
 * Stores mesh in 2 .obj files:
 *   robin_fuselage
 *   robin_pylon
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
  // when c[7] (C8 in the paper) is zero, this will throw NaNs
  return c[5] + c[6]*std::pow(c[0]+c[1]*std::pow((x+c[2])/c[3], c[4]), 1.0/c[7]);
}

double getRadialCoord(double H, double W, double theta, double N) {
  double numer = 0.25*H*W;
  // Without std::abs() here, this line throws NaNs
  double denom = std::pow(0.5*H*std::sin(theta), N) + std::pow(0.5*W*std::cos(theta), N);
  return numer / std::pow(denom, 1.0/N); 
}

double chebeshev_node(double a, double b, double k, double n) {
  if (k == 0) {
    return a;
  } else if (k==n) {
    return b;
  } else {
    return (a+b)*0.5 + (b-a)*0.5*std::cos((2.0*(n-k))*M_PI*0.5/n);
  }
}

void create_vertices(const size_t nx, const size_t nt, const std::string fileName,
                     int (*getSection)(double), const double xBegin, const double xEnd) {

  // first 4 rows of each section are the fuselage
  // next 2 more rows each for the pylon
  // These are somewhat corrected coefficients from 1987
  std::vector<SupEll> hcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, 0.0, 0.25, 1.8},
                                 {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
                                 {1.0, -1.0, -0.8, 0.4, 3.0, 0.0, 0.145, 3.0},
                                 {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.145, 2.0} };

  std::vector<SupEll> wcoeff = { {1.0, -1.0, -0.4, 0.4, 2.0, 0.0, 0.25, 2.0},
                                 {0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.05, 0.2, 0.6},
                                 {1.0, -1.0, -1.9, 0.1, 2.0, 0.0, 0.05, 2.0},
                                 {1.0, -1.0, -0.8, 0.4, 3.0, 0.0, 0.166, 3.0},
                                 {1.0, -1.0, -0.8, 0.218, 2.0, 0.0, 0.166, 2.0} };

  std::vector<SupEll> zcoeff = { {1.0, -1.0, -0.4, 0.4, 1.8, -0.08, 0.08, 1.8},
                                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.04, -0.04, 0.6},
                                 {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {1.0, -1.0, -0.8, 1.1, 1.5, 0.065, 0.06, 0.6} };

  std::vector<SupEll> ncoeff = { {2.0, 3.0, 0.0, 0.4, 1.0, 0.0, 1.0, 1.0},
                                 {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {5.0, -3.0, -0.8, 1.1, 1.0, 0.0, 0.0, 0.0},
                                 {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                 {5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };

  // Open file to write to
  std::ofstream file;
  file.open(fileName);
  file << "# Vertices\n";

  std::cout << "  generating nodes" << std::endl;
  for (size_t ix=0; ix<nx+1; ix++) {
    
    const double xol = chebeshev_node(xBegin, xEnd, ix, nx);
    //std::cout << ix << " " << xol << std::endl;
    const int sec = getSection(xol);
    if (sec == -1) {
      std::cout << "ERROR: sec == " << sec << " x == " << xol << std::endl;
      exit(0);
    }

    // compute H, W, Z0, and N from xol and the constants
    const double H  = getsuperval(xol, hcoeff[sec]);
    const double W  = getsuperval(xol, wcoeff[sec]);
    const double Z0 = getsuperval(xol, zcoeff[sec]);
    const double N  = getsuperval(xol, ncoeff[sec]);

    for (size_t it=0; it<nt; it++) {
      const double theta = 2.0*M_PI*it/(double)nt;
      // compute r from H, W, N, theta
      const double r = getRadialCoord(H, W, theta, N);
      // compute yol, zol from r, theta, Z0
      double yol = r * std::sin(theta);
      double zol = r * std::cos(theta) + Z0;
      //std::cout << xol << " " << yol << " " << zol << std::endl;

      // change nans to 0.0 so something shows up
      if (std::isnan(yol)) yol = 0.0;
      if (std::isnan(zol)) zol = 0.0;
      
      // Write vertex to file
      file << "v " << xol << " " << yol << " " << zol << "\n";
      if ((ix == 0) || (ix == nx)) { break; }
    }
  }
  file.close();
}

void create_faces(const std::string fileName, const size_t nx, const size_t nt) {
  std::ofstream file;
  file.open(fileName, std::ios_base::app);
  std::cout << "  generating faces" << std::endl;

  // Label faces
  file << "\n# Faces\n";
  for (size_t i=2; i<nt+1; i++) {
    file << "f " << 1 << " " << i << " " << i+1 << "\n";
  }
  file << "f " << 1 << " " << nt+1 << " " << 2 << "\n";

  // Link the two rings together to make faces
  for (size_t r=0; r<nx-2; r++) {
    // node indices for rings 1 and 2
    const size_t ir1 = nt*r + 2;
    const size_t ir2 = nt*(r+1) + 2;

    // Runs through all nodes in the ring
    for (size_t n=0; n<nt-1; n++) {
      file << "f " << ir1+n << " " << ir2+n << " " << ir1+n+1 << "\n";
      file << "f " << ir2+n << " " << ir2+n+1 << " " << ir1+n+1 << "\n";
    }
    file << "f " << ir1+nt-1 << " " << ir2+nt-1 << " " << ir1 << "\n";
    file << "f " << ir2+nt-1 << " " << ir2 << " " << ir1 << "\n";
  }

  size_t lastVert = nt*(nx-1)+2;
  for (size_t i=2; i<nt+1; i++) {
    file << "f " << (i+1)+nt*(nx-2) << " " << i+nt*(nx-2) << " " << lastVert << "\n";
  }
  file << "f " << 2+nt*(nx-2) << " " << lastVert-1 << " " << lastVert << "\n";
  // Done writing faces
  file.close();
}

long long input_check(std::string num) {
  long long n;
  try {
    std::size_t pos;
    n = std::stoll(num, &pos);
    // check if input is mix of numbers and chars
    // 1zz converts to 1
    if (pos < num.size()) {
      std::cerr << "Trailing characters after number: " << num << '\n';
      std::cerr << "Using: " << n << std::endl;
    } 
  } catch (std::invalid_argument const &ex) {
    std::cerr << "Invalid number: " << num << '\n';
    exit(0);
  } catch (std::out_of_range const &ex) {
    std::cerr << "Number out of range: " << num << '\n';
    exit(0);
  }
  return n;
}

// execution starts here
int main(int argc, char const *argv[]) {
  std::cout << "Generating ROBIN model" << std::endl;

  // Actual coefficients are in create_vertices

  // Read in command line inputs
  if (argc != 5) {
    std::cerr << "Incorrect number of arguments: " << argc-1 << std::endl;
    std::cerr << "Argument order is [nx fuselage] [nt fuselage] [nx pylon] [nt pylon]\n";
    exit(0);
  }
  const long long fnx = input_check((std::string)argv[1]);
  const long long fnt = input_check((std::string)argv[2]);
  const long long pnx = input_check((std::string)argv[3]);
  const long long pnt = input_check((std::string)argv[4]);
 
  const std::string fileName = "robin_fuselage.obj"; 
  const double fusBegin = 0.0;
  const double fusEnd = 2.0;
  const std::string pFileName = "robin_pylon.obj";
  const double pylBegin = 0.4;
  const double pylEnd = 1.018;

  std::cout << "Creating fuselage mesh..." << std::endl;
  create_vertices(fnx, fnt, fileName, get_fuselage_section, fusBegin, fusEnd);
  create_faces(fileName, fnx, fnt);
  std::cout << "Creating pylon mesh..." << std::endl;
  create_vertices(pnx, pnt, pFileName, get_pylon_section, pylBegin, pylEnd);
  create_faces(pFileName, pnx, pnt);
  // CSG these two meshes together 
  return 0;
}
