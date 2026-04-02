//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sn_passive.cpp
//! \brief Initializes SN ejecta into polar wedge

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

// outflow condition with diode
void SNOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// inject SN ejecta
void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// set time step
Real MyTimeStep(MeshBlock* pmb);

namespace {
// problem parameters which are useful to make global to this file
Real rho0, p0, gammagas, gasEntropy;
bool diode;
Real Rej, Eej, Mej, vmax;
Real smallDt, largeDt, smallDtDuration;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

// read parameters from input file
void Mesh::InitUserMeshData(ParameterInput *pin) {

  // ambient medium properties
  rho0       = pin->GetOrAddReal("problem","rho0",1.0);
  p0         = pin->GetOrAddReal("problem","p0",1.0);

  // information for thermal energy
  gammagas   = pin->GetReal("hydro","gamma");
  gasEntropy = pin->GetReal("problem","gasEntropy");

  // use diode condition for outflow
  diode      = pin->GetOrAddBoolean("problem","diode",false);

  // properties of ejecta
  Rej        = pin->GetOrAddReal("problem","Rejecta",1.0);
  Eej        = pin->GetOrAddReal("problem","Eejecta",1.0);
  Mej        = pin->GetOrAddReal("problem","Mejecta",1.0);
  vmax       = pin->GetOrAddReal("problem","vmax",1.0);

  // properties of initial timestep
  smallDt = pin->GetReal("problem", "smallDt");
  largeDt = pin->GetReal("problem", "largeDt");
  smallDtDuration = pin->GetReal("problem", "smallDtDuration");

  // tell athena to use our user-defined boundaries
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);

  // set a very small timestep initially
  EnrollUserTimeStepFunction(MyTimeStep);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes SN ejecta.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

Real r, t0, v0sq, v, rho, pres, prefactor;

// Extract floors
Real rho_floor = peos->GetDensityFloor();
Real p_floor = peos->GetPressureFloor();

t0 = Rej/vmax;
v0sq = 4.0/3.0*Eej/Mej;
prefactor = std::pow(3.0/4.0/3.14159/Eej, 1.5) * std::pow(Mej, 2.5);

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);

        // transverse velocities are always zero
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (r > Rej || true) {
          // ambient medium
          phydro->u(IDN,k,j,i) = rho0;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = p0/(gammagas-1.0); // + 0.5*rho0*vel0*vel0;
        } else {
          // ejecta profile
          v = vmax*r/Rej;
          rho = prefactor * std::exp( -v*v/v0sq ) /t0/t0/t0;
          pres = gasEntropy * std::pow(rho,1.6666666666667);
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = rho*v;
          phydro->u(IEN,k,j,i) = pres/(gammagas-1.0) + 0.5*rho*v*v;

        //initialize everything at the beginning to be 1.0 for the passive tracers, to get all the original material
        pscalars->r(0, k, j, i) = 1.0;
        pscalars->r(1, k, j, i) = 0.0;

        }

      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SNOuterX1()
//  \brief Sets boundary condition on downstream outflow boundary

void SNOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
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
//  \brief Sets inner boundary where SN ejecta is injected

void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rho, vel, v0sq, pres, prefactor, t0, t_fix;
  // Extract floors
  Real rho_floor = pmb->peos->GetDensityFloor();
  Real p_floor = pmb->peos->GetPressureFloor();

  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  t0 = Rej/vmax;
    //choose constant t instead of using time and t0
    //be the Gaussian at one point in time, what number to choose, time that gives velocity of ~10k or ~20k 
  t_fix = Rej * ((1.0/(1.5*std::pow(10, 9))) - (1.0/vmax)); //solve vel equation for time for 20k km/s
  //t_fix = Rej * ((1.0/(8.5*std::pow(10, 8))) - (1.0/vmax)); //8.5k km/s
  //t_fix = Rej * ((1.0/(1.2*std::pow(10, 8))) - (1.0/vmax)); //12k km/s
  v0sq = 4.0/3.0*Eej/Mej;

  vel = 1.0/(t_fix/Rej + 1.0/vmax);
  prefactor = std::pow(3.0/4.0/3.14159/Eej,1.5)*std::pow(Mej,2.5);
  rho = prefactor*std::exp(-vel*vel/v0sq)*std::pow(t_fix+t0,-3.0);
  pres = gasEntropy * std::pow(rho,1.6666666666667);

  if ( Rej==0.0 || Eej==0.0 || Mej==0.0 || vmax==0.0 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in sn.cpp ProblemGenerator" << std::endl
        << "ejecta parameters not set" << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = rho;
        prim(IM1,k,j,il-i) = vel;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pres;


         //if time is less than 1 crossing times keep the first passive tracer one and the second zero  
        if (time < 15.0) { //15 for 20k km/s, 35 for 8.5k km/s, 25 for 12k
            prim_scalar(0,k,j,il-i) = 1.0;
            prim_scalar(1,k,j,il-i) = 0.0;
        }
            //otherwise make the first one zero and the second one 1 
        else {
            prim_scalar(0,k,j,il-i) = 0.0;
            prim_scalar(1,k,j,il-i) = 1.0;

        }


     }
    }
  }
}

Real MyTimeStep(MeshBlock *pmb) {
  Real time   = pmb->pmy_mesh->time;
  Real dt     = pmb->pmy_mesh->dt;
  Real min_dt = largeDt;

  if ( time < smallDtDuration ) {
    min_dt = std::min( dt, smallDt );
  }

  return min_dt;
}
