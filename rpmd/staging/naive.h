#ifndef naive_h
#define naive_h

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

    random_device rd;    //random number device
    mt19937 gen;         //random number generator, need to set up seed in init()


//--------------- FUNCTIONS -----------------
Real delta(int a, int b);

Real heaviside(Real x);

Real get_wall_time();

Real get_cpu_time();


//---------------- PROPAGATION ----------------
void INIT_HAMILTONIAN();
void MOVE_VV_V(Real dt);
void MOVE_VV_R(Real dt);
void UPDATE_RP_F();
void UPDATE_EX_F();
void MOVE_LF_V(Real dt);
void MOVE_LF_R(Real dt);
void INIT_SAMPLING_Wigner();
void INIT_SAMPLING_Cl();
void PropOneTrajVV(int t);
void PropOneTrajLF();
void STG2NAIVE();    //Naive to Staging variables
void NAIVE2STG();    //Staging Variables
void STGRPF();    //Update Staging force
void STGRPP(Real dt); //Update RP momenta
void STGRPR(Real dt); //Update position
void MOVE_STGV(Real dt); //Update External influence on velocity

//---------------- Variables ----------------
vector<Real> R;    //nuclear position coordinate
vector<Real> V;    //nuclear velocity
vector<Real> P;    //nuclear momenta
vector<Real> F;    //nuclear force
vector<Real> F1;   //Internal Ring Polymer (RP) force
vector<Real> stgR; //Staging Position Coordinates
vector<Real> stgV; //Staging Velocity Coordinates
vector<Real> stgP; //Staging Momenta 
vector<Real> stgF; //Staging Force (inner RP)

vector<Real> mass;    //nuclear mass
vector<Real> vmass;
vector<Real> ovmass;  //one over mass
vector<Real> ovvmass;

Real Energy;
Real KEne;
Real PEne;
vector<Real> E;
vector<Real> U;
vector<Real> KE;
vector<Real> U0;
vector<Real> U1;
vector<Real> U2;
vector<Real> U3;
vector<Real> KE0;
vector<Real> KE1;
vector<Real> KE2;
vector<Real> KE3;

vector<int> position;
vector<int> displacement;
vector<int> velocity;


//---------------- CONSTANTS ----------------
Real DT = 0.00001;
Real DT2 = DT * 0.5;
const Real hbar = 1;
const Real pi = std::acos(-1.0);//3.14159265358979;
const Complex I(0.0,1.0);  //the imaginary I
int DOFn = 1;          // How many nuclear Degree of Freedom(DOFn) are talking about
int nsteps = 10000000;        // Total MD propagation Steps
int sampling = 100;
int nbeads = 16;            // RPMD beads number
Real beta = 1.0;
Real temp = 1.0 / beta;
vector<Real> Req;      // Equilibrated coordinates
vector<Real> Omega;    // Force field frequency
vector<Real> w2;       // Omega^2
const Real wp = sqrt(static_cast<Real> (nbeads))* temp / hbar; //RP Quadro const
const Real wp2 = wp * wp;

#endif
