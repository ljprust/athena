//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file octree_test.cpp
//! \brief Problem generator for reading in data using an octree search

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

// need C environment so Fortran subroutine names aren't garbled
extern "C" {

namespace{
  extern void build_tree(int *numRead, Real *boxSize);
  extern void tree_walk(Real *x, Real *y, Real *z,
                        int *id, Real *boxSize);
  extern void octree_final();
}

int NumToRead;
Real boxSize, gammagas;

void Mesh::InitUserMeshData(ParameterInput *pin) {

  NumToRead = pin->GetInteger("problem","NumToRead");
  boxSize = pin->GetReal("mesh","x1max");
  gammagas = pin->GetReal("hydro","gamma");

  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  build_tree(&NumToRead, &boxSize);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem generator for testing the octree search
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  std::vector<Real> rho_in, vr_in, pres_in;
  Real rho, vr, pres;
  Real r, theta, x, y, z;
  int index;
  char rhoFile[256];
  char vrFile[256];
  char presFile[256];

  sprintf(rhoFile,"rho.txt");
  sprintf(vrFile,"vr.txt");
  sprintf(presFile,"pres.txt");

  std::ifstream rhoFileRead, vrFileRead, presFileRead;
  rhoFileRead.open(rhoFile);
  vrFileRead.open(vrFile);
  presFileRead.open(presFile);

  for (int l=0; l<NumToRead; l++) {
    rhoFileRead  >> rho;
    vrFileRead   >> vr;
    presFileRead >> pres;

    rho_in.push_back(rho);
    vr_in.push_back(vr);
    pres_in.push_back(pres);
  }
  rhoFileRead.close();
  vrFileRead.close();
  presFileRead.close();

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        theta = pcoord->x2v(i);
        x = r*std::cos(theta);
        y = r*std::sin(theta);
        z = 0.0;

        tree_walk(&x, &y, &z, &index, &boxSize);

        phydro->u(IDN,k,j,i) = rho_in[index];
        phydro->u(IM1,k,j,i) = vr_in[index];
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres_in[index]/(gammagas-1.0)
          + 0.5*rho_in[index]*vr_in[index]*vr_in[index];
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Deallocate octree
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  octree_final();
  return;
}

} // extern C
