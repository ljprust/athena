//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file mesa.cpp
//! \brief implements ideal EOS in general EOS framework, mostly for debuging
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../eos.hpp"

extern "C" {

namespace{
  const Real a_rad = 7.5646e-15; // ergs / (cm^3 K^4)
  const Real R_gas = 8.314e7; // ergs / (mol K)
  bool debug = true;
  char MesaDir[256] = "/Users/ljprust/code/mesa-r10398";

  // set the composition
  Real X = 1.0;
  Real Z = 0.0;
  int use_solar = 1;
  Real fc12 = 0.0; Real fn14 = 0.0; Real fo16 = 0.0; Real fne20 = 0.0;

  extern void mesaeos_DTget(Real *Rho, Real *T, Real *Xin,
    Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
    Real *press, Real *energy, Real *gamma);

  extern void mesaeos_DEget(Real *Rho, Real *energy, Real *T_guess, Real *Xin,
    Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
    Real *T, Real *press, Real *gamma);

  // extern void mesaeos_get_T_given_Ptotal(Real *Rho, Real *T_guess, Real *Xin,
  //   Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16,
  //   Real *fne20, Real *T, Real *press, Real *gamma);

  extern void mesaeos_init(char MesaDir);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return total pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {

  // we will get these from DEget
  Real press, gamma, T, pres;

  // guess T assuming ideal monatomic gas
  Real Tguess = egas/(3.0/2.0*rho*R_gas);

  // call MESA
  mesaeos_DEget( &rho, &egas, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &pres, &gamma);

  return pres;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {

  // we will get these from MESA
  Real egas, gamma, T;

  // guess E and T assuming ideal monatomic gas
  Real Eguess = 3.0/2.0*pres;
  Real Tguess = Eguess/(3.0/2.0*rho*R_gas);

  // call MESA
  mesaeos_DEget( &rho, &Eguess, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &pres, &gamma);

  // call MESA
  //mesaeos_get_T_given_Ptotal( &rho, &T, &X, &Z, &use_solar, &fc12,
  //  &fn14, &fo16, &fne20, &pres, &egas, &gamma);

  // plug T back in to get E, gamma
  mesaeos_DTget( &rho, &T, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &pres, &egas, &gamma);

  return egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {

  // we will get these from MESA
  Real egas, gamma, T;

  // guess E assuming ideal monatomic gas
  Real Tguess = pres/(rho*R_gas);
  Real Eguess = 3.0/2.0*pres;

  // call MESA
  mesaeos_DEget( &rho, &Eguess, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &pres, &gamma);

  // call MESA
  //mesaeos_get_T_given_Ptotal( &rho, &T, &X, &Z, &use_solar, &fc12,
  //  &fn14, &fo16, &fne20, &pres, &egas, &gamma);

  // plug T back in to get E, gamma
  mesaeos_DTget( &rho, &T, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &pres, &egas, &gamma);

  return gamma * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  mesaeos_init( *MesaDir );
  return;
}
}
