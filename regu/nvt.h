#ifndef nvt_h
#define nvt_h

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
void UPDATE_F();
void MOVE_OU(Real dt);
void MOVE_LF_V(Real dt);
void MOVE_LF_R(Real dt);
void INIT_SAMPLING_Wigner();
void INIT_SAMPLING_Cl();
void PropOneTrajVV(int t);
void PropOneTrajLF();
void INIT_Langevin();
void INIT_NHC();
void MOVE_NHC(Real dt);
void MOVE();
void PrepOneTrajVV(int t);
void MOVE_Nose(Real dt);
void INIT_Nose();
void MOVE_Nose_V(Real dt);
void MOVE_Nose_R(Real dt);
void MOVE_NoseOU(Real dt);
void INIT_Regu();
void MOVE_ReguR(Real dt);
void MOVE_ReguP(Real dt);
void MOVE_Reguv1(Real dt);
void MOVE_ReguOU(Real dt);
void MOVE_ReguSC(Real dt);

//---------------- Variables ----------------
vector<Real> R;    //nuclear position coordinate
vector<Real> V;    //nuclear velocity 
vector<Real> mass;    //nuclear mass
vector<Real> ovmass;    //one over mass
vector<Real> P;    //nuclear momenta
vector<Real> F;    //nuclear force

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
vector<Real> Q;       //Nose Hoover thermostat mass
vector<Real> ovQ;     // 1/Q
vector<Real> fric;    //Langevin thermostat friction term 
vector<Real> sigma;   //Ornstein-Ohlenbeck process distribution STD
vector<Real> v;    //Thermostat momenta
vector<Real> G;    //Thermostat "Force"
vector<Real> wsy;  //Suzuki-Yolanda Factorisation
vector<Real> tsy;  //Suzuki-Yolanda Factorisation 
vector<Real> v1;    // SIN-R 1st thermostat momenta
vector<Real> Q1;    // SIN-R 1st thermostat masses
vector<Real> ovQ1;  // 1/Q1
vector<Real> v2;    // SIN-R 2nd thermostat momenta
vector<Real> Q2;    // SIN-R 2nd thermostat masses
vector<Real> ovQ2;  // 1/Q2
vector<Real> c;
vector<Real> ovc;


long int position[10000];
long int displacement[10000];
long int velocity[10000];


Real DT;    //molecular dynamic time interval
Real DT2;   //half of MD time interval
Real beta;  //beta = 1/ k_BT, inverse temperature
Real temp;  //Temperature T 
//---------------- CONSTANTS ----------------
const Real hbar = 1;
const Real pi = std::acos(-1.0);//3.14159265358979;
const Complex I(0.0,1.0);  //the imaginary I
int DOFn = 2;          // How many nuclear Degree of Freedom(DOFn) are talking about
int nsteppr = 100000;       // MD prep MD steps 
int nsteps = 1000000;        // Propagation MD propagation Steps
int sampling = 1;
int nres = 5;          // RESPA loop number
int nsy = 3;           // Forth order Suzuki-Yolanda factorisation
int NHL = 10;
vector<Real> Req;      // Equilibrated coordinates
vector<Real> Omega;    // Force field frequency

#endif
