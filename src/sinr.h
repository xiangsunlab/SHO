#ifndef sinr_h
#define sinr_h

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

//-----------------Step-related parameters------------------------
const int nstep = 10;  //MD propagation steps
const int nres  = 5;   //RESPA factorisation are done in nres steps
const int nsy = 3;     //4th order Suzuki-Yolanda factorisatio number are done in 3 steps
const int NHL = 4;     //There are NHL sets of theromstat extension co-ordinates or momenta

//-----------------MD related constants---------------------------
const real hbar = 1.0;  //Reduced Planck constant
const real k = 1.0;     //Boltzmann constant
const real M = 1.0;     //Mass for the particle
const real beta = 1.0;  //Inversed temperature
const real kBT = 1.0/beta; //Temperature*Boltzmann Constant for the system 

//-----------------MD trajectory related variables---------------
real R;    //Position coordinate
real V;    //Velocity V=dR/dt; Consider V as 1st-order time derivation for the coordinate
real P;    //Momenta  P=dH/d xdot
real F;    //Force    F=-dH/dR
real E;    //System total energy
real KE;   //System total kinetic energy
real PE;    //System total potential

//-----------------Thermostat Constant---------------------------
const real fric = 0.1;    //Friction term in the Langevin thermostat
const real Q1 = 1.0;      //Pseudo-Mass for the Nose-Hoover-like thermostat for the first coupling layer
const real Q2 = 1.0;      //Pseudo-Mass for the Nose-Hoover-like thermostat for the second coupling layer
const real sigma = sqrt(2.0 * fric / (beta * Q2));  //Sigma for the Ornstein-Uhlenbeck process

//-----------------Function for MD---------------------------------
void addup(real dt);
#endif
