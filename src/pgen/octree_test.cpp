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
std::vector<Real> rho_in, vr_in, pres_in;

void Mesh::InitUserMeshData(ParameterInput *pin) {

  NumToRead = pin->GetInteger("problem","NumToRead");
  boxSize = pin->GetReal("mesh","x1max");
  gammagas = pin->GetReal("hydro","gamma");
/*
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }
*/
  char rhoFile[256];
  char vrFile[256];
  char presFile[256];

  sprintf(rhoFile,"rho.txt");
  sprintf(vrFile,"vr.txt");
  sprintf(presFile,"pres.txt");

  printf("opening data files\n");
  std::ifstream rhoFileRead, vrFileRead, presFileRead;
  rhoFileRead.open(rhoFile);
  vrFileRead.open(vrFile);
  presFileRead.open(presFile);

  Real rho, vr, pres;
  for (int l=0; l<NumToRead; l++) {
    rhoFileRead  >> rho;
    vrFileRead   >> vr;
    presFileRead >> pres;

    rho_in.push_back(rho);
    vr_in.push_back(vr);
    pres_in.push_back(pres);
  }
  printf("done reading, closing data files\n");
  rhoFileRead.close();
  vrFileRead.close();
  presFileRead.close();
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem generator for testing the octree search
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real r, theta, x, y, z;
  int index;

  printf("Building tree\n");
  build_tree(&NumToRead, &boxSize);

  printf("beginning treewalk\n");
  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        theta = pcoord->x2v(j);
        x = r; // r*std::cos(theta);
        y = theta; // r*std::sin(theta);
        z = 0.0;

        index = 0;
        tree_walk(&x, &y, &z, &index, &boxSize);
        // if (i==is && j==js) index = 1; // hack for float min bug
        index = index-1; // different indexing in C++ vs Fortran
        printf("for r theta %5.3e %5.3e found neighbor id %d with rho %5.3e\n",r,theta,index,rho_in[index]);

        phydro->u(IDN,k,j,i) = rho_in[index];
        phydro->u(IM1,k,j,i) = vr_in[index];
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres_in[index]/(gammagas-1.0)
          + 0.5*rho_in[index]*vr_in[index]*vr_in[index];
      }
    }
  }
  printf("done with treewalk\n");
  //printf("Calling octree_final\n");
  //octree_final();
  //printf("done with octree_final\n");
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Deallocate octree
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  printf("Calling octree_final\n");
  octree_final();
  return;
}

} // extern C
