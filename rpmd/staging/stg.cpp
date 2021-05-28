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
void MOVE_VV_V(Real dt);  // Naive Propagation of Velocity under Internal And EXTERNAL FF
void MOVE_VV_R(Real dt);  // Naive Propagation of Position after velocity update
void UPDATE_RP_F();       // Naive Update of the External FF
void NAIVE2STG();         // Cartisian coordinates to staging variables
void STG2NAIVE();         // Staging coordinates to Cartisian coordinates
void UPDATE_STG_RPF();    // Staging Update of Ring Polymer Force
void UPDATE_NM_RPF();     // Normal Mode Update of Ring Polymer Force
void MOVE_STG_RPV(Real dt);      // Staging Velocity Propagation of Ring Polymer 
void MOVE_NM_RPV(Real dt);       // Normal Mode Velocity Propagation of Ring Polymer
void MOVE_STG_RPR(Real dt);      // Staging Position Propagation of Ring Polymer Force
void MOVE_NM_RPR(Real dt);       // Normal Mode Position Propagation of Ring Polymer
void UPDATE_EX_F();       // Naive Propagation of External FF
void MOVE_STG_EXF();      // Staging Propagation of External FF   
void MOVE_LF_V(Real dt);
void MOVE_LF_R(Real dt);
void INIT_SAMPLING_Wigner();
void INIT_SAMPLING_Cl();
void PropOneTrajVV(int t);
void PropOneTrajLF();

//---------------- Variables ----------------
vector<Real> R;    //nuclear position Cartisian coordinate
vector<Real> RP_R; //nuclear position STAGING coordinate
vector<Real> V;    //nuclear Cartisian velocity
vector<Real> RP_V; //nuclear STAGING velocity
vector<Real> mass;    //nuclear Real mass
vector<Real> vmass;   //Staging Mass
vector<Real> ovmass;  //one over Real mass
vector<Real> ovvmass; //one over Staging Mass
vector<Real> P;       //nuclear Cartisian momenta
vector<Real> RP_P;    //nuclear Ring Polymer momenta
vector<Real> F;    //nuclear Cartisian force from external Force Field(FF)
vector<Real> F1;   //Ring Polymer force

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
int DOFn = 2;          // How many nuclear Degree of Freedom(DOFn) are talking about
int nsteps = 1;        // Total MD propagation Steps
int sampling = 1;
int nbeads = 16;            // RPMD beads number
Real beta = 1.0;
Real temp = 1.0 / beta;
vector<Real> Req;      // Equilibrated coordinates
vector<Real> Omega;    // Force field frequency
vector<Real> w2;       // Omega^2
const Real wp = sqrt(static_cast<Real> (nbeads))* temp / hbar; //RP Quadro const
const Real wp2 = wp * wp;

#endif
