//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ideal.cpp
//! \brief implements ideal EOS in general EOS framework, mostly for debuging
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../eos.hpp"

Real pressureEq(Real temp, Real rho);
Real energyEq(Real temp, Real rho);
Real calcGamma(Real rho, Real pres);
Real TempFromBeta3Gamma4(Real beta3, Real gamma4);

namespace{
  const Real a       = 7.5646e-15;
  const Real mu      = 2.0;
  const Real mProton = 1.6726e-24;
  const Real kB      = 1.3807e-16;
}

//takes in temperature and density and outputs pressure
Real pressureEq(Real temp, Real rho) {
  return 1.0/3.0*a*std::pow(temp,4.0) + rho*kB*temp/(mu*mProton);
}

//takes in temperature and density and outputs internal energy
Real energyEq(Real temp, Real rho) {
  return 3.0/2.0*kB*temp*rho/mu/mProton + a*temp*temp*temp*temp;
}

//takes in temperature pressure and density and returns gamma
Real calcGamma(Real rho, Real pres) {
  Real gasPres, beta;

  Real beta3 = 3.0*rho*kB/mu/mProton/a;
  Real gamma4 = 3.0*pres/a;
  Real temp = TempFromBeta3Gamma4(beta3, gamma4);

  gasPres = rho*kB*temp/(mu*mProton);
  beta = gasPres/pres;
  return (32.0-24.0*beta-3.0*beta*beta)/(24.0-21.0*beta);
}

Real TempFromBeta3Gamma4(Real beta3, Real gamma4) {
  Real temp;
  Real beta = std::pow(beta3, 1.0/3.0);
  Real gamma = std::pow(gamma4, 0.25);
  Real epsilon = beta/gamma;
  Real epsilon3 = beta*beta*beta/gamma/gamma/gamma;
  Real delta = 1.0/epsilon;
  if( epsilon < 1.0 ) {
    temp = gamma*( 1.0 - 0.25*epsilon3 
         - 1.0/32.0*epsilon3*epsilon3
         + 7.0/2048.0*epsilon3*epsilon3*epsilon3*epsilon3
         + 1.0/512.0*epsilon3*epsilon3*epsilon3*epsilon3*epsilon3 );
  } else {
    temp = beta*( std::pow(delta,4.0) - std::pow(delta,16.0)
                + 4.0*std::pow(delta,28.0) - 22.0*std::pow(delta,40.0) );
  }
  return temp;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  Real beta3 = 3.0/2.0*rho*kB/mu/mProton/a;
  Real gamma4 = egas/a;
  Real temp = TempFromBeta3Gamma4(beta3, gamma4);
  //printf("PresFromRhoEg: %5.3e\n", temp);
  return pressureEq(temp, rho);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  Real beta3 = 3.0*rho*kB/mu/mProton/a;
  Real gamma4 = 3.0*pres/a;
  Real temp = TempFromBeta3Gamma4(beta3, gamma4);
  //printf("EgasFromRhoP: %5.3e\n", temp);
  return energyEq(temp, rho);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  Real gamma1=calcGamma(rho,pres);
  //printf("gamma = %5.3e\n",gamma1);
  return gamma1 * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  return;
}
