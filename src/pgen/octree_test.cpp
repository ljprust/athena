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
  extern void build_tree(int *numRead, Real *x1max);
  extern void tree_walk(Real *x, Real *y, // Real *z,
                        int *id, Real *deltax);
  extern void octree_final();
}

int NumToRead, nx1;
Real x1min, x1max, deltax, gammagas;
std::vector<Real> rho_in, vr_in, pres_in;

void Mesh::InitUserMeshData(ParameterInput *pin) {

  NumToRead = pin->GetInteger("problem","NumToRead");
  x1min = pin->GetReal("mesh","x1min");
  x1max = pin->GetReal("mesh","x1max");
  nx1 = pin->GetInteger("mesh","nx1");
  gammagas = pin->GetReal("hydro","gamma");

  deltax = (x1max-x1min)/(static_cast<double>(nx1)); // average particle spacing

  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in octree_test.cpp ProblemGenerator" << std::endl
        << "Spherical-polar coordainates are assumed: " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  char rhoFile[256];
  char vrFile[256];
  char presFile[256];

  sprintf(rhoFile,"rho.txt");
  sprintf(vrFile,"vr.txt");
  sprintf(presFile,"pres.txt");

  printf("Opening data files with state variables...\n");
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
  printf("Done reading, closing data files\n");
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

  Real r, theta, x, y;
  int index;

  printf("Building tree...\n");
  build_tree(&NumToRead, &x1max);
  printf("Octree build success!\n");

  printf("Beginning tree walk...\n");
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        theta = pcoord->x2v(j);
        // ignoring phi on purpose
        x = r*std::cos(theta);
        y = r*std::sin(theta);
        //printf("x y boxsize %5.3e %5.3e %5.3e\n",x,y,x1max);

        index = 0;
        tree_walk(&x, &y, &index, &deltax);
        // if (i==is && j==js) index = 1; // hack for double min bug
        index = index-1; // correct for C++ vs Fortran indices
        // printf("for r theta %5.3e %5.3e found neighbor id %d with rho %5.3e\n",r,theta,index,rho_in[index]);

        phydro->u(IDN,k,j,i) = rho_in[index];
        phydro->u(IM1,k,j,i) = rho_in[index]*vr_in[index];
        phydro->u(IM2,k,j,i) = 0.0; // assuming homology so v_theta = v_phi = 0
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres_in[index]/(gammagas-1.0)
          + 0.5*rho_in[index]*vr_in[index]*vr_in[index];
        // do we need to call the EOS here?
      }
    }
  }
  printf("Finished tree walk\n");
  printf("Dealloating octree...\n");
  octree_final();
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Deallocate octree
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  //printf("Calling octree_final\n");
  //octree_final();
  return;
}

} // extern C
