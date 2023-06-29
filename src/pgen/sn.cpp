//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sn.cpp
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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// inflow/outflow BCs
void SNOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// problem parameters which are useful to make global to this file
Real rho0, p0, gammagas;
bool diode;
Real Rej, Eej, Mej, vmax;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  rho0     = pin->GetOrAddReal("problem","rho0",1.0);
  p0       = pin->GetOrAddReal("problem","p0",1.0);
  gammagas = pin->GetOrAddReal("hydro","gamma",0.0);
  diode    = pin->GetOrAddBoolean("problem","diode",false);
  Rej      = pin->GetOrAddReal("problem","Rejecta",1.0);
  Eej      = pin->GetOrAddReal("problem","Eejecta",1.0);
  Mej      = pin->GetOrAddReal("problem","Mejecta",1.0);
  vmax     = pin->GetOrAddReal("problem","vmax",1.0);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
  //EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes SN ejecta.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

Real r, t0, v0, v, rho, pres;

t0 = Rej/vmax;
v0 = std::sqrt(4.0/3.0*Eej/Mej);

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);

        // transverse velocities are always zero
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (r > Rej) {
          // ambient medium
          phydro->u(IDN,k,j,i) = rho0;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = p0/(gammagas-1.0); // + 0.5*rho0*vel0*vel0;
        } else {
          // ejecta profile
          v = vmax*r/Rej;
          rho = std::pow(3.0/4.0/3.14159/Eej,1.5)*std::pow(Mej,2.5)
                *std::exp(-v*v/v0/v0)/t0/t0/t0;
          pres = std::pow(rho,gammagas);
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = v;
          phydro->u(IEN,k,j,i) = pres/(gammagas-1.0) + 0.5*rho*v*v;
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
/*
  Real rho, pres;

  if ( pvacuum==0.0 || dvacuum==0.0 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
        << "vacuum pressure and/or density not set" << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = dvacuum;
        prim(IM1,k,j,il-i) = 0.0;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pvacuum;

      }
    }
  }
*/
}

//namespace {}
