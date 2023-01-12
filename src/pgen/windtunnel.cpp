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

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
// problem parameters which are useful to make global to this file
Real gm0, rho0, gamma_gas;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel, p0;
  Real x1, x2, x3;
  Real gm1 = peos->GetGamma() - 1.0;

  den = 1.0;
  vel = 1.0;
  p0  = 1.0;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        phydro->u(IDN,k,j,i) = den;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM1,k,j,i) =  vel*std::cos(phi); // radial
          phydro->u(IM2,k,j,i) = -vel*std::sin(phi); // azimuth
          phydro->u(IM3,k,j,i) =  0.0;               // z
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM1,k,j,i) =  vel*std::cos(phi); // radial
          phydro->u(IM2,k,j,i) = -vel*std::sin(phi); // polar
          phydro->u(IM3,k,j,i) =  0.0;               // azimuth
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          phydro->u(IM1,k,j,i) =  vel; // x
          phydro->u(IM2,k,j,i) =  0.0; // y
          phydro->u(IM3,k,j,i) =  0.0; // z
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "no acceptable coord system found" << std::endl;
          ATHENA_ERROR(msg);
        }

        /*
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

        phydro->u(IEN,k,j,i) = p0/gm1
          + (0.5)*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
          + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        */

        phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*vel*vel / den;
      }
    }
  }

  return;
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
