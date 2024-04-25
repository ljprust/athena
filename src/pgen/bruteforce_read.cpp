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

int NumToRead;
Real gammagas, x1max;
std::vector<Real> r_in, theta_in, x_in, y_in, rho_in, vr_in, pres_in;

void Mesh::InitUserMeshData(ParameterInput *pin) {

  NumToRead = pin->GetInteger("problem","NumToRead");
  gammagas = pin->GetReal("hydro","gamma");
  x1max = pin->GetReal("mesh","x1max");

  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in bruteforce_read.cpp ProblemGenerator" << std::endl
        << "Spherical-polar coordainates are assumed: " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  char rFile[256];
  char thetaFile[256];
  char rhoFile[256];
  char vrFile[256];
  char presFile[256];

  sprintf(rFile,    "r.txt");
  sprintf(thetaFile,"theta.txt");
  sprintf(rhoFile,  "rho.txt");
  sprintf(vrFile,   "vr.txt");
  sprintf(presFile, "pres.txt");

  printf("Opening data files with state variables...\n");
  std::ifstream rFileRead, thetaFileRead, rhoFileRead, vrFileRead, presFileRead;
  rFileRead.open(rFile);
  thetaFileRead.open(thetaFile);
  rhoFileRead.open(rhoFile);
  vrFileRead.open(vrFile);
  presFileRead.open(presFile);

  Real r, theta, rho, vr, pres;
  for (int l=0; l<NumToRead; l++) {
    rFileRead     >> r;
    thetaFileRead >> theta;
    rhoFileRead   >> rho;
    vrFileRead    >> vr;
    presFileRead  >> pres;

    r_in.push_back(r);
    theta_in.push_back(theta);
    x_in.push_back( r * std::cos(theta) );
    y_in.push_back( r * std::sin(theta) );
    rho_in.push_back(rho);
    vr_in.push_back(vr);
    pres_in.push_back(pres);
  }
  printf("Done reading, closing data files\n");
  rFileRead.close();
  thetaFileRead.close();
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

  Real r, theta, x, y, dist2, minDist2;
  int index;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        theta = pcoord->x2v(j);
        // ignoring phi on purpose
        x = r*std::cos(theta);
        y = r*std::sin(theta);

        index = -1;
        minDist2 = x1max*x1max;
        for (int l=0; l<NumToRead; l++) {
          dist2 = (x-x_in[l])*(x-x_in[l]) + (y-y_in[l])*(y-y_in[l]);
          if ( dist2 < minDist2 ) {
            minDist2 = dist2;
            index = l;
          }
        }
        //printf("for r theta %5.3e %5.3e found neighbor id %d with rho %5.3e\n",r,theta,index,rho_in[index]);

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
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Deallocate octree
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}
