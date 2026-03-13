//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file windtunnel.cpp
//! \brief Initializes parallel flow in one direction in both cylindrical and
//! spherical polar coordinates.

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
void BinaryWind2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                         Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void BinaryWind2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                         Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
// problem parameters which are useful to make global to this file
//Real gm_companion, r_companion; 
Real separation, gammagas, Mdot_wind, v_wind, pressure_ratio, r_inner;
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
  r_inner = pin->GetOrAddReal("mesh","x1min",0.0);
  //r_companion = pin->GetOrAddReal("problem","r_companion",0.0);
  separation = pin->GetOrAddReal("problem","separation",0.0);
  gammagas = pin->GetOrAddReal("hydro","gamma",0.0);
  Mdot_wind = pin->GetOrAddReal("problem","Mdot_wind",0.0);
  v_wind = pin->GetOrAddReal("problem","v_wind",0.0); 
  pressure_ratio = pin->GetOrAddReal("problem","pressure_ratio",0.0);
  diode = pin->GetOrAddBoolean("problem","diode",false);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, BinaryWind2DOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, BinaryWind2DInnerX1);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes wind tunnel.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rho, x1;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    //x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      //x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
/*
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          y = x1*std::sin(x2);
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          y = x1*std::sin(x2)*std::cos(x3);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          y = x2;
        }
*/
        rho = Mdot_wind/v_wind/4.0/3.14159/x1/x1;

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
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          // THIS IS NOT SUPPORTED YET!!!
          phydro->u(IM1,k,j,i) = rho*v_wind; // x
          phydro->u(IM2,k,j,i) = 0.0; // y
          phydro->u(IM3,k,j,i) = 0.0; // z
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

void BinaryWind2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;
  //Real phi;
  //Real rho, pres;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        //rad=pco->x1v(iu+i);
        //phi=pco->x2v(j);
        //z=pco->x3v(k);
/*
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          y = pco->x2v(iu+i)*std::sin(pco->x1v(j));
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          y = pco->x1v(iu+i)*std::sin(pco->x2v(j))*std::cos(pco->x3v(k));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          y = pco->x2v(j);
        }
*/
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

void BinaryWind2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rho, pres;

  rho = Mdot_wind/4.0/3.14159/v_wind/r_inner/r_inner;
  pres = pressure_ratio*0.5*rho*v_wind*v_wind;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = rho;
        prim(IM1,k,j,il-i) = v_wind;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pres;

      }
    }
  }
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

}
