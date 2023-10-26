//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shock_tube.cpp
//! \brief Problem generator for shock tube problems.
//!
//! Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
//! shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), freopen(), fprintf(), fclose()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

extern "C" {

namespace {
Real rho0, vel0, pres0, gammagas;
Real a_rad;
Real R_gas;

// set the composition
//Real X;
//Real Z;
//int use_solar; 
//Real fc12, fn14, fo16 = 0.0, fne20 = 0.0;

extern void mesaeos_dtget(Real *Rho, Real *T, Real *Xin,
  Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
  Real *press, Real *energy, Real *gamma);

extern void mesaeos_deget(Real *Rho, Real *energy, Real *T_guess, Real *Xin,
  Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
  Real *T, Real *press, Real *gamma);

extern void mesaeos_dtget_t_given_ptotal(Real *Rho, Real *T_guess, Real *press,
  Real *Xin, Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16,
  Real *fne20, Real *T, Real *gamma);

extern void mesaeos_init();

} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  rho0  = pin->GetOrAddReal("problem","rho0",1.0);
  pres0 = pin->GetOrAddReal("problem","pres0",1.0);
  vel0  = pin->GetOrAddReal("problem","vel0",0.0);
  gammagas = pin->GetOrAddReal("hydro","gamma",0.0);

  printf("Calling mesa eos init from pgen ...\n");
  mesaeos_init();
/*
  a_rad = 7.5646e-15; // ergs / (cm^3 K^4)
  R_gas = 8.314e7; // ergs / (mol K)

  // set the composition
  Real X = 1.0;
  Real Z = 0.0;
  int use_solar = 1;  
  Real fc12 = 0.0; Real fn14 = 0.0; Real fo16 = 0.0; Real fne20 = 0.0;
*/
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  a_rad = 7.5646e-15; // ergs / (cm^3 K^4)
  R_gas = 8.314e7; // ergs / (mol K)

  // set the composition
  Real X = 1.0;
  Real Z = 0.0;
  int use_solar = 1;
  Real fc12 = 0.0; Real fn14 = 0.0; Real fo16 = 0.0; Real fne20 = 0.0;

  Real T, gamma, presJunk, Especific;

  Real Tguess = std::min(std::pow(3.0*pres0/a_rad, 0.25), pres0/R_gas/rho0);
  printf("Guessing T = %5.3e\n",Tguess);

  mesaeos_dtget_t_given_ptotal( &rho0, &Tguess, &pres0, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &gamma);
  printf("intermediate T gamma from rho Tguess pres %5.3e %5.3e %5.3e %5.3e %5.3e\n",T,gamma,rho0,Tguess,pres0);
  mesaeos_dtget( &rho0, &T, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);
  printf("initial energy density = %5.3e\n",Especific*rho0);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = rho0;
        phydro->u(IM1,k,j,i) = rho0*vel0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        //phydro->u(IEN,k,j,i) = pres0/(gammagas-1.0) + 0.5*rho0*vel0*vel0;
        phydro->u(IEN,k,j,i) = Especific*rho0 + 0.5*rho0*vel0*vel0;
      }
    }
  }
  return;
}
}
