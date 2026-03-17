//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binarywind.cpp
//! \brief Initializes a wind from a central star within the inner radial boundary.

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
void BinaryWindOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void BinaryWindInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
//Real gm_companion, r_companion; 
Real Mdot_wind, v_wind, pressure_ratio;
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
  //gm_primary = pin->GetOrAddReal("problem","GM",0.0);
  //gm_companion = pin->GetOrAddReal("problem","gm_companion",0.0);
  //r_companion = pin->GetOrAddReal("problem","r_companion",0.0);
  Mdot_wind = pin->GetOrAddReal("problem","Mdot_wind",0.0);
  v_wind = pin->GetOrAddReal("problem","v_wind",0.0); 
  pressure_ratio = pin->GetOrAddReal("problem","pressure_ratio",0.0);
  diode = pin->GetOrAddBoolean("problem","diode",false);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, BinaryWindOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, BinaryWindInnerX1);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes wind tunnel.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rho, r;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    //x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      //x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);

        rho = Mdot_wind/v_wind/4.0/3.14159/r/r;

        phydro->u(IDN,k,j,i) = rho;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          //GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
          phydro->u(IM1,k,j,i) = rho*v_wind; // radial
          phydro->u(IM2,k,j,i) = 0.0;        // azimuth
          phydro->u(IM3,k,j,i) = 0.0;        // z
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM1,k,j,i) = rho*v_wind; // radial
          phydro->u(IM2,k,j,i) = 0.0;        // polar
          phydro->u(IM3,k,j,i) = 0.0;        // azimuth
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "no acceptable coord system found" << std::endl;
          ATHENA_ERROR(msg);
        }

        phydro->u(IEN,k,j,i) = pressure_ratio*0.5*rho*v_wind*v_wind;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DOuterX1()
//  \brief Sets boundary condition on upstream boundary (oib) for wind tunnel
//
// Quantities at this boundary are held fixed at the constant upstream state

void BinaryWindOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

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
//! \fn void WindTunnel2DInnerX1()
//  \brief Sets vacuum inner boundary
//
// Quantities in ghost cells are set to some pressure and density

void BinaryWindInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real r, rho, pres;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        r = pco->x1v(il-i);
        rho = Mdot_wind/4.0/3.14159/v_wind/r/r;
        pres = pressure_ratio*0.5*rho*v_wind*v_wind;

        prim(IDN,k,j,il-i) = rho;
        prim(IM1,k,j,il-i) = v_wind;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pres;

      }
    }
  }
}

//namespace {
//}
