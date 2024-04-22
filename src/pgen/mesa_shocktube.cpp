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
Real rho0L, vel0L, pres0L;
Real rho0R, vel0R, pres0R;
Real gammagas, a_rad, R_gas;

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
  rho0L    = pin->GetOrAddReal("problem","rho0L", 1.0);
  pres0L   = pin->GetOrAddReal("problem","pres0L",1.0);
  vel0L    = pin->GetOrAddReal("problem","vel0L", 0.0);
  rho0R    = pin->GetOrAddReal("problem","rho0R", 1.0);
  pres0R   = pin->GetOrAddReal("problem","pres0R",1.0);
  vel0R    = pin->GetOrAddReal("problem","vel0R", 0.0);
  gammagas = pin->GetOrAddReal("hydro",  "gamma", 0.0);

  printf("Calling mesa eos init from pgen ...\n");
  mesaeos_init();

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  a_rad = 7.5646e-15; // ergs / (cm^3 K^4)
  R_gas = 8.314e7; // ergs / (mol K)
  Real X = 1.0;
  Real Z = 0.0;
  int use_solar = 1;
  Real fc12 = 0.0; Real fn14 = 0.0; Real fo16 = 0.0; Real fne20 = 0.0;
  Real x;

  Real TL, gammaL, presJunkL, EspecificL;
  Real TR, gammaR, presJunkR, EspecificR;

  Real TguessL = std::min(std::pow(3.0*pres0L/a_rad, 0.25), pres0L/R_gas/rho0L);
  printf("Guessing TL = %5.3e\n",TguessL);
  mesaeos_dtget_t_given_ptotal( &rho0L, &TguessL, &pres0L, &X, &Z,
    &use_solar, &fc12, &fn14, &fo16, &fne20, &TL, &gammaL);
  printf("intermediate T gamma from rhoL TguessL presL %5.3e %5.3e %5.3e %5.3e %5.3e\n",
  TL,gammaL,rho0L,TguessL,pres0L);
  mesaeos_dtget( &rho0L, &TL, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunkL, &EspecificL, &gammaL);

  Real TguessR = std::min(std::pow(3.0*pres0R/a_rad, 0.25), pres0R/R_gas/rho0R);
  printf("Guessing TR = %5.3e\n",TguessR);
  mesaeos_dtget_t_given_ptotal( &rho0R, &TguessR, &pres0R, &X, &Z,
    &use_solar, &fc12, &fn14, &fo16, &fne20, &TR, &gammaR);
  printf("intermediate T gamma from rhoR TguessR presR %5.3e %5.3e %5.3e %5.3e %5.3e\n",
  TR,gammaR,rho0R,TguessR,pres0R);
  mesaeos_dtget( &rho0R, &TR, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunkR, &EspecificR, &gammaR);

  //EspecificL = pres0L/rho0L/(5.0/3.0-1.0);
  //EspecificR = pres0R/rho0R/(5.0/3.0-1.0);
  printf("initial energy density left  = %5.3e\n",EspecificL*rho0L);
  printf("initial energy density right = %5.3e\n",EspecificR*rho0R);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        x = pcoord->x1v(i);
        if (x<0.0) {
          phydro->u(IDN,k,j,i) = rho0L;
          phydro->u(IM1,k,j,i) = rho0L*vel0L;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          //phydro->u(IEN,k,j,i) = pres0/(gammagas-1.0) + 0.5*rho0*vel0*vel0;
          phydro->u(IEN,k,j,i) = EspecificL*rho0L + 0.5*rho0L*vel0L*vel0L;
        } else {
          phydro->u(IDN,k,j,i) = rho0R;
          phydro->u(IM1,k,j,i) = rho0R*vel0R;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = EspecificR*rho0R + 0.5*rho0R*vel0R*vel0R;
        }
      }
    }
  }
  return;
}
}
