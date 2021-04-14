#ifndef aux_h
#define aux_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <random>
#include <complex>
#include <iomanip>
#include <stdexcept>
#include <memory>
#include <sys/time.h>
using namespace std;


//--------------- TYPES -------------------
typedef double Real;
typedef std::complex<Real>  Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<vector<Complex> > > Complex_3D_Array;
typedef std::vector<vector<Real> > Real_Matrix;
typedef std::vector<vector<vector<Real> > > Real_3D_Array;


//--------------- FUNCTIONS -----------------
Real delta(int a, int b);

Real heaviside(Real x);

Real get_wall_time();

Real get_cpu_time();




//---------------- CONSTANTS ----------------
const Real hbar = 1;
const Real pi = std::acos(-1.0);//3.14159265358979;
const Complex I(0.0,1.0);  //the imaginary I


