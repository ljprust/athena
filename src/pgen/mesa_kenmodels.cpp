//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shock_tube.cpp
//! \brief Problem generator for shock tube problems.
//!
//! Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
//! shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), freopen(), fprintf(), fclose()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

extern "C" {

namespace {
Real Mej, Eej, r0, rhoFloor, vmax;
Real X, Z, fc12, fn14, fo16, fne20;
int use_solar;
int NumToRead;
std::vector<Real> rho_in, vr_in, temp_in;

extern void mesaeos_dtget(Real *Rho, Real *T, Real *Xin,
  Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
  Real *press, Real *energy, Real *gamma);

extern void mesaeos_deget(Real *Rho, Real *energy, Real *T_guess, Real *Xin,
  Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16, Real *fne20,
  Real *T, Real *press, Real *gamma);

extern void mesaeos_dtget_t_given_ptotal(Real *Rho, Real *T_guess, Real *press,
  Real *Xin, Real *Zin, int *use_solar, Real *fc12, Real *fn14, Real *fo16,
  Real *fne20, Real *T, Real *gamma);

extern void mesaeos_init(const char *MesaDir);

} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Mej       = pin->GetReal("problem","Mej");
  Eej       = pin->GetReal("problem","Eej");
  r0        = pin->GetReal("problem","r0");
  rhoFloor  = pin->GetReal("problem","rhoFloor");
  Tfloor    = pin->GetReal("problem","Tfloor");
  NumToRead = pin->GetInteger("problem","NumToRead");

  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in bruteforce_read.cpp ProblemGenerator" << std::endl
        << "Spherical-polar coordainates are assumed: " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  char rhoFile[256];
  char vrFile[256];
  char tempFile[256];

  sprintf(rhoFile,  "rho.txt");
  sprintf(vrFile,   "vr.txt");
  sprintf(tempFile, "temp.txt");

  printf("Opening data files with state variables...\n");
  std::ifstream rhoFileRead, vrFileRead, tempFileRead;
  rhoFileRead.open(rhoFile);
  vrFileRead.open(vrFile);
  tempFileRead.open(tempFile);

  Real rho, vr, temp;
  vmax = 0.0;
  for (int l=0; l<NumToRead; l++) {
    rhoFileRead   >> rho;
    vrFileRead    >> vr;
    tempFileRead  >> temp;

    if( vr > vmax ) vmax = vr;

    rho_in.push_back(rho);
    vr_in.push_back(vr);
    temp_in.push_back(temp);
  }
  printf("Done reading, closing data files\n");
  rhoFileRead.close();
  vrFileRead.close();
  tempFileRead.close();

  std::string MesaDir = pin->GetString("problem","MesaDir");
  char cMesaDir[256];
  strncpy(cMesaDir, MesaDir.c_str(), sizeof(cMesaDir));
  cMesaDir[sizeof(cMesaDir) - 1] = 0;

  X = pin->GetOrAddReal("problem","X",1.0);
  Z = pin->GetOrAddReal("problem","Z",0.0);
  use_solar = pin->GetOrAddInteger("problem","use_solar",1);

  // metal compositions as a fraction of Z
  fc12  = pin->GetOrAddReal("problem","fc12", 0.5);
  fn14  = pin->GetOrAddReal("problem","fn14", 0.0);
  fo16  = pin->GetOrAddReal("problem","fo16", 0.5);
  fne20 = pin->GetOrAddReal("problem","fne20",0.0);

  printf("\nCalling MESA EOS init from problem generator ...\n");
  mesaeos_init(cMesaDir);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real r, t0, vr, rho;
  Real temp, gamma, presJunk;
  Real Especific, EspecificFloor;
  Real dv, mindv;
  int index;

  t0 = r0/vmax;

  mesaeos_dtget( &rhoFloor, &Tfloor, &X, &Z, &use_solar, &fc12,
    &fn14, &fo16, &fne20, &presJunk, &EspecificFloor, &gamma);

  printf("initial energy density floor = %5.3e\n",EspecificFloor*rhoFloor);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);
        if ( x < r0 ) {
          vr = r/t0;

          index = -1;
          mindv = 10.0e9;
          for (int l=0; l<NumToRead; l++) {
            dv = abs(vr-vr_in[l]);
            if ( dv < mindv ) {
              mindv = dv;
              index = l;
            }
          }
          //printf("for r %5.3e found neighbor id %d with rho %5.3e\n",r,index,rho_in[index]);

          rho = rho_in[index];
          vr = vr_in[index];
          temp = temp_in[index];

          mesaeos_dtget( &rho, &temp, &X, &Z, &use_solar, &fc12,
            &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);

          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = rho*vr;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = Especific*rho + 0.5*rho*vr*vr;
        } else {
          phydro->u(IDN,k,j,i) = rhoFloor;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = EspecificFloor*rhoFloor; // + 0.5*rho0R*vel0R*vel0R;
        }
      }
    }
  }
  return;
}
}
