//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cee_sn.cpp
//! \brief Initializes CEE outflows within domain and sets SN inner radial boundary.

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// inflow/outflow BCs
void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

Real MyTimeStep(MeshBlock* pmb);

void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt, 
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar);

namespace {
Real gammagas, Rgas, vmax, ramPressureFactor, rhoISM, r_inner, Rsun, Mej, Eej, t0;
Real Mdotwind, vwind;
Real day;
Real epsilon_Ni, epsilon_Co;
Real tau_Ni, tau_Co;
Real A_nuc, mproton;
Real index_Ni, index_Co;
bool diode;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gammagas          = pin->GetOrAddReal("hydro","gamma",0.0);
  diode             = pin->GetOrAddBoolean("problem","diode",false);
  vmax              = pin->GetOrAddReal("problem","v_max",0.0);
  ramPressureFactor = pin->GetOrAddReal("problem","ramPressureFactor",0.0);
  rhoISM            = pin->GetOrAddReal("problem","rho_ISM",0.0);
  r_inner           = pin->GetOrAddReal("mesh","x1min",0.0);
  Mej               = pin->GetOrAddReal("problem","Mej",0.0);
  Eej               = pin->GetOrAddReal("problem","Eej",0.0);
  Rsun              = 7.0e10;
  Rgas              = 8.314e7;
  Mdotwind          = 1.0e-8*2.0e33/365.25/24.0/3600.0;
  vwind             = 30.0e5;
  t0                = r_inner/vmax;
  day               = 24.0*3600.0;
  mproton           = 1.6726e-24;
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiodeOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);
  EnrollUserTimeStepFunction(MyTimeStep);
  //EnrollUserExplicitSourceFunction(radioactiveHeating);

  epsilon_Ni = 1.72e-6; // energy released in decays
  epsilon_Co = 3.49e-6;
  tau_Ni     = 8.77*day; // e-folding times of decays
  tau_Co     = 111.0*day;
  A_nuc      = 56.0; // mass number

  index_Ni = 0;
  index_Co = 0;

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes CEE outflows within domain.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real r, theta, z;
  Real diskHeight, rhoCEE, rhoMin, rhoWind;
  Real rho, temp, pres, mintemp;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    //phi = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      theta = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);
        z = r*std::cos(theta);

        diskHeight = 0.4*r+10.0*Rsun;
		// Rsun*(95.0*std::log10(r/Rsun)-125.0);
	rhoCEE = std::exp(-z*z/2.0/diskHeight/diskHeight)
	//* 0.01*std::pow(r/10.0/Rsun,-4.0)*std::pow(1.0+std::pow(125.0*Rsun/r,3.5),-1.05); // 500 d
	* (1.0e-4*std::pow(r/100.0/Rsun,-4.8)*std::pow(1.0+std::pow(690.0*Rsun/r,1.9),-3.4)+6.358e-13*std::pow(r/6000.0/Rsun,-3.2)/(1.0+std::pow(6000.0*Rsun/r,6.0))); // 10000 d
	
	//rhoMin = 1.0e-11;
	rhoWind = Mdotwind/4.0/3.14159/r/r/vwind;

	mintemp = 1.0e4;

	if (rhoCEE > rhoWind && rhoCEE > rhoISM) { // disk
	  rho = rhoCEE;
    temp = std::max( mintemp, 4.5e4/(r/100.0/Rsun) );
    pres = rho*Rgas*temp;
	} else if (rhoWind > rhoISM) { // wind
	  rho = rhoWind;
    pres = rho*Rgas*mintemp;
	} else { // ISM
    rho = rhoISM;
    pres = ramPressureFactor*rhoISM*vmax*vmax;
  }

        phydro->u(IDN,k,j,i) = rho;

        phydro->u(IM1,k,j,i) = 0.0; // rho*vel0*std::cos(x2); // radial
        phydro->u(IM2,k,j,i) = 0.0; //-rho*vel0*std::sin(x2); // polar
        phydro->u(IM3,k,j,i) = 0.0;               // azimuth

        pscalars->s(0,k,j,i) = 0.0;

        phydro->u(IEN,k,j,i) = pres;
		            // ramPressureFactor*rhoCEE*vmax*vmax; 
                            // pres/(gammagas-1.0) + 0.5*rho*vel0*vel0;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void DiodeOuterX1()
//  \brief Sets doide outflow conditions at outer x1 boundary
//
// Gas outflows with optional diode condition

void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        // let gas outflow freely
        prim(IDN,k,j,iu+i) = prim(IDN,k,j,iu);
        prim(IM2,k,j,iu+i) = prim(IM2,k,j,iu);
        prim(IM3,k,j,iu+i) = prim(IM3,k,j,iu);
        prim(IEN,k,j,iu+i) = prim(IEN,k,j,iu);

        // ensure that no gas enters through the outflow boundary
        // by giving the radial velocity some TLC
        applyDiode = prim(IM1,k,j,iu) < 0.0;
        if (diode && applyDiode) {
          prim(IM1,k,j,iu+i) = 0.0;
        } else {
          prim(IM1,k,j,iu+i) = prim(IM1,k,j,iu);
        }

      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void SNInnerX1()
//  \brief Sets supernova ejecta inner boundary
//
// Quantities in ghost cells are set to Gaussian ejecta model from Wong+24

void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rhoSunny, pres;
  Real v0sq, v_inner, prefactor;

  v_inner = r_inner/(time + t0);
  v0sq = 4.0 / 3.0 * Eej / Mej;

  prefactor = std::pow(3.0 / 4.0 / 3.14159 / Eej, 1.5) * std::pow(Mej, 2.5);
  rhoSunny = prefactor * std::exp(-v_inner * v_inner / v0sq) * std::pow(time + t0, -3.0);
  pres = 0.7e14*std::pow(rhoSunny,1.6666666666667);
  //pres = ramPressureFactor*rhoSunny*vmax*vmax;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i)        = rhoSunny;
        prim(IM1,k,j,il-i)        = v_inner;
        prim(IM2,k,j,il-i)        = 0.0;
        prim(IM3,k,j,il-i)        = 0.0;
        prim(IEN,k,j,il-i)        = pres;
        pmb->pscalars->r(0,k,j,i) = 1.0;

      }
    }
  }
}

Real MyTimeStep(MeshBlock* pmb) {
    Real time   = pmb->pmy_mesh->time;
    Real dt     = pmb->pmy_mesh->dt;
    Real min_dt = 1.0e3;
    if (time < 0.1) min_dt = std::min(dt, 1.0e-2);
    return min_dt;
}

void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt,
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar) {

    Real deltaRho_Ni, deltaRho_Co;

    Real coefficient_Ni = epsilon_Ni / (tau_Ni * A_nuc * mproton);
    Real coefficient_Co = epsilon_Co / (tau_Co * A_nuc * mproton);

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
            for (int i = pmb->is; i <= pmb->ie; ++i) {
                cons(IEN,k,j,i) += dt * cons_scalar(index_Ni,k,j,i) * coefficient_Ni 
                                 * std::exp(-1.0 * (time + t0) / tau_Ni)
                                 + dt * cons_scalar(index_Co,k,j,i) * coefficient_Co
                                 * std::exp(-1.0 * (time + t0) / tau_Co);
                deltaRho_Ni = -dt/tau_Ni*cons_scalar(index_Ni,k,j,i);
                deltaRho_Co = -dt/tau_Co*cons_scalar(index_Co,k,j,i);
                cons_scalar(index_Ni,k,j,i) += deltaRho_Ni;
                cons_scalar(index_Co,k,j,i) += deltaRho_Co - deltaRho_Ni;
            }
        }
    }
}  

namespace {
}
