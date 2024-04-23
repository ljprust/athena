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
Real vmax, Mej, Eej, r0, rhoFloor, presFloor;
Real a_rad, R_gas;

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
  vmax      = pin->GetReal("problem","vmax");
  Mej       = pin->GetReal("problem","Mej");
  Eej       = pin->GetReal("problem","Eej");
  r0        = pin->GetReal("problem","r0");
  rhoFloor  = pin->GetReal("problem","rhoFloor");
  //presFloor = pin->GetReal("problem","presFloor");

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
  Real x, t0, v0sq, prefactor, vr, rho, pres;

  t0 = r0/vmax;
  v0sq = 4.0/3.0*Eej/Mej;
  prefactor = std::pow(3.0/4.0/3.14159/Eej,1.5)*std::pow(Mej,2.5)/t0/t0/t0;

  Real temp, Tguess, Tfloor, gamma, presJunk, Especific, EspecificFloor;

  //Tguess = std::min(std::pow(3.0*presFloor/a_rad, 0.25), presFloor/R_gas/rhoFloor);
  //printf("Guessing T = %5.3e\n",Tguess);
  //mesaeos_dtget_t_given_ptotal( &rhoFloor, &Tguess, &presFloor, &X, &Z,
  //  &use_solar, &fc12, &fn14, &fo16, &fne20, &temp, &gamma);
  //printf("intermediate T gamma from rhoL TguessL presL %5.3e %5.3e %5.3e %5.3e %5.3e\n",
  //TL,gammaL,rho0L,TguessL,pres0L);
  Tfloor = 1.0e5;
  mesaeos_dtget( &rhoFloor, &Tfloor, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunk, &EspecificFloor, &gamma);

  //EspecificL = pres0L/rho0L/(5.0/3.0-1.0);
  printf("initial energy density floor = %5.3e\n",EspecificFloor*rhoFloor);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        x = pcoord->x1v(i);
        if (x<r0) {
          vr = x/t0;
          rho = prefactor*std::exp(-vr*vr/v0sq);

          temp = 1.0e6; // at t0 = 2 hr
          //Tguess = std::min(std::pow(3.0*presFloor/a_rad, 0.25), presFloor/R_gas/rhoFloor);
          //printf("Guessing T = %5.3e\n",Tguess);
          //mesaeos_dtget_t_given_ptotal( &rhoFloor, &Tguess, &presFloor, &X, &Z,
          //  &use_solar, &fc12, &fn14, &fo16, &fne20, &temp, &gamma);
          //printf("intermediate T gamma from rhoL TguessL presL %5.3e %5.3e %5.3e %5.3e %5.3e\n",
          //TL,gammaL,rho0L,TguessL,pres0L);
          mesaeos_dtget( &rho, &temp, &X, &Z, &use_solar, &fc12,
            &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);

          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = rho*vr;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          //phydro->u(IEN,k,j,i) = pres0/(gammagas-1.0) + 0.5*rho0*vel0*vel0;
          phydro->u(IEN,k,j,i) = Especific*rho + 0.5*rho*vr*vr;
        } else {
          phydro->u(IDN,k,j,i) = rhoFloor;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = EspecificFloor*rhoFloor; // + 0.5*rho0R*vel0R*vel0R;
        }
      }
    }
  }
  return;
}
}
