//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file mesa.cpp
//! \brief implements MESA EOS
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../eos.hpp"

extern "C" {

namespace{
  const Real a_rad = 7.5646e-15; // ergs / (cm^3 K^4)
  const Real R_gas = 8.314e7; // ergs / (mol K)
  bool debug = false;
  char MesaDir[256] = "/Users/ljprust/code/mesa-r10398";

  // set the composition
  Real X = 1.0;
  Real Z = 0.0;
  int use_solar = 1;
  Real fc12 = 0.0; Real fn14 = 0.0; Real fo16 = 0.0; Real fne20 = 0.0;

  // __mesaeos_lib_MOD_mesaeos_dtget( ...
  extern void mesaeos_dtget(Real *Rho, Real *T, Real *Xin,
    Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
    Real *press, Real *energy, Real *gamma);

  extern void mesaeos_deget(Real *Rho, Real *energy, Real *T_guess, Real *Xin,
    Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
    Real *T, Real *press, Real *gamma);

  extern void mesaeos_dtget_t_given_ptotal(Real *Rho, Real *T_guess, Real *press,
    Real *Xin, Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16,
    Real *fne20, Real *T, Real *gamma);

  //extern void mesaeos_init(char MesaDir);
  extern void mesaeos_init();
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return total pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {

  // we will get these from DEget
  Real gamma, T, pres;

  // convert to specific energy
  Real Especific = egas/rho;

  // guess temperature assuming it's either
  // radiation- or gas-pressure dominated
  Real Tguess = std::min(Especific*(5.0/3.0-1.0)/R_gas, std::pow(egas/a_rad, 0.25));

  // MESA EOS call
  mesaeos_deget( &rho, &Especific, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &pres, &gamma);

  if (debug) printf("egas->pres: rho egas Tguess T pres gamma %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n",rho,egas,Tguess,T,pres,gamma);

  return pres;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {

  // we will get these from MESA
  Real egas, Especific, gamma, T;

  // we won't need the pressure we get from the EOS call
  Real presJunk;

  //if(rho<0.0) printf("NEGATIVE DENSITY IN D,P -> E!!! \n");

  // guess temperature assuming it's either
  // radiation- or gas-pressure dominated
  Real Tguess = std::min(std::pow(3.0*pres/a_rad, 0.25), pres/R_gas/rho);

  // first call EOS to get temperature
  mesaeos_dtget_t_given_ptotal( &rho, &Tguess, &pres, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &gamma);

  if(debug) printf("intermediate T gamma from rho Tguess pres %5.3e %5.3e %5.3e %5.3e %5.3e\n",T,gamma,rho,Tguess,pres);

  // then use the temperature to get specific energy
  mesaeos_dtget( &rho, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);

  if(debug) printf("pres->egas: rho pres Tguess presJunk egas gamma %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n",rho,pres,Tguess,presJunk,Especific*rho,gamma);

  // convert back to energy density
  return Especific*rho;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {

  // we will get these from MESA
  Real egas, Especific, gamma, T, cs2;

  // we won't need the pressure we get from the EOS call
  Real presJunk;

  // guess temperature assuming it's either
  // radiation- or gas-pressure dominated
  Real Tguess = std::min(std::pow(3.0*pres/a_rad, 0.25), pres/R_gas/rho);

  // first call EOS to get temperature
  mesaeos_dtget_t_given_ptotal( &rho, &Tguess, &pres, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &T, &gamma);

  if(debug) printf("intermediate T gamma from rho Tguess pres %5.3e %5.3e %5.3e %5.3e %5.3e\n",T,gamma,rho,Tguess,pres);

  // then use temperature to get gamma
  mesaeos_dtget( &rho, &Tguess, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);
  // we ignore the pressure and energy outputs from this call

  if(debug) printf("for cs2, using gamma = %5.3e\n",gamma);

  cs2 = gamma*pres/rho;

  return cs2;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  //mesaeos_init( *MesaDir );
  //mesaeos_init();
  MesaDir = GetString("problem","MesaDir");

  X = GetOrAddReal("problem","X",1.0);
  Z = GetOrAddReal("problem","Z",0.0);

  // flag to use solar metal fractions
  use_solar = GetOrAddInteger("problem","use_solar",1);

  // metal compositions as a fraction of Z
  fc12  = GetOrAddReal("problem","fc12", 0.5);
  fn14  = GetOrAddReal("problem","fn14", 0.0);
  fo16  = GetOrAddReal("problem","fo16", 0.5);
  fne20 = GetOrAddReal("problem","fne20",0.0);

  return;
}
} // extern C
