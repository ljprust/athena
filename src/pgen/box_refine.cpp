//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file box.cpp
//! \brief Initializes boring uniform box

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

Real BrokenGeometric(Real x, RegionSize rs);

namespace {
Real rho0, vel0, p0, gammagas;
int Ntot, Nin;
Real Rbreak, ratio;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  rho0     = pin->GetOrAddReal("problem","rho0",1.0);
  vel0     = pin->GetOrAddReal("problem","vel0",1.0);
  p0       = pin->GetOrAddReal("problem","p0",1.0);
  gammagas = pin->GetOrAddReal("hydro","gamma",0.0);

  Rbreak = pin->GetOrAddReal("mesh","Rbreak",0.0);
  ratio  = pin->GetOrAddReal("mesh","ratio",0.0);
  Nin    = pin->GetOrAddInteger("mesh","Nin",0);
  Ntot   = pin->GetOrAddInteger("mesh","nx1",0);
  //EnrollUserMeshGenerator(X1DIR, BrokenGeometric);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes wind tunnel.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = rho0;
        phydro->u(IM1,k,j,i) = rho0*vel0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = p0/(gammagas-1.0) + 0.5*rho0*vel0*vel0;

        pscalars->s(0,k,j,i) = 0.0;
      }
    }
  }

  return;
}

Real BrokenGeometric(Real x, RegionSize rs) {
  printf("here 0\n");
  int index = static_cast<int>(x * static_cast<double>(Ntot))-1;
  printf("here 1\n");
  std::vector<Real> xout;
  std::vector<Real> deltaxout;
  Real deltaxin = (Rbreak-rs.x1min)/static_cast<double>(Nin);
  printf("here 2\n");
  Real ratiosum = 0.0;
  for (int i=Nin; i<Ntot; ++i) { ratiosum += std::pow(ratio, i-Nin); }
  printf("here 3\n");
  Real deltaxbreak = (rs.x1max-Rbreak)/ratiosum;
  for (int j=0; j<Ntot; ++j) {
    if(j<Nin) {
      deltaxout.push_back(deltaxin);
      xout.push_back(rs.x1min + deltaxin*static_cast<double>(j) + 0.5*deltaxin);
    } else {
      deltaxout.push_back(deltaxbreak*std::pow(ratio, j-Nin));
      xout.push_back(xout[j-1] + 0.5*deltaxout[j-1] + 0.5*deltaxout[j]);
    }
  }
  printf("here 4\n");
  return xout[index];
}
