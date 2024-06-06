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
#include <fstream>
#include <iostream>

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
Real r0, rhoFloor, Tfloor, vmax, tdata, gammaGas;
int NumToRead;
std::vector<Real> rho_in, vr_in, temp_in;
bool useGammaLaw;

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
  r0        = pin->GetReal("problem","r0");
  rhoFloor  = pin->GetReal("problem","rhoFloor");
  Tfloor    = pin->GetReal("problem","Tfloor");
  NumToRead = pin->GetInteger("problem","NumToRead");
  useGammaLaw = pin->GetOrAddBoolean("problem","useGammaLaw");
  gammaGas  = pin->GetOrAddReal("hydro","gamma",1.666666667);

  // time of data in input files, needed for rescaling density
  // negative if no difference from t0
  tdata     = pin->GetOrAddReal("problem","tdata",-1.0);

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

  if (!useGammaLaw) {
    printf("\nCalling MESA EOS init from problem generator ...\n");
    mesaeos_init(cMesaDir);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real r, t0, vr, rho;
  //Real xc, yc, zc;
  Real temp, gamma, presJunk;
  Real Especific, EspecificFloor;
  Real dv, mindv, scaleFactorRho, scaleFactorTemp;
  int index;

  Real X, Z, fc12, fn14, fo16, fne20;
  int use_solar;

  X = pin->GetOrAddReal("problem","X",1.0);
  Z = pin->GetOrAddReal("problem","Z",0.0);
  use_solar = pin->GetOrAddInteger("problem","use_solar",1);

  // metal compositions as a fraction of Z
  fc12  = pin->GetOrAddReal("problem","fc12", 0.5);
  fn14  = pin->GetOrAddReal("problem","fn14", 0.0);
  fo16  = pin->GetOrAddReal("problem","fo16", 0.5);
  fne20 = pin->GetOrAddReal("problem","fne20",0.0);

  // center of SN
  //xc = 0.5*pin->GetReal("mesh","x1max");
  //yc = 0.5*pin->GetReal("mesh","x2max");
  //zc = 0.5*pin->GetReal("mesh","x3max");

  t0 = r0/vmax;

  if ( tdata > 0.0 ) {
    scaleFactorRho  = std::pow( tdata/t0, 3.0 );
    scaleFactorTemp = std::pow( tdata/t0, 2.0 );
  } else {
    scaleFactorRho  = 1.0;
    scaleFactorTemp = 1.0;
  }

  if(!useGammaLaw) {
    mesaeos_dtget( &rhoFloor, &Tfloor, &X, &Z, &use_solar, &fc12,
      &fn14, &fo16, &fne20, &presJunk, &EspecificFloor, &gamma);
    printf("initial energy density floor = %5.3e\n",EspecificFloor*rhoFloor);
  }

  for (int k=ks; k<=ke; ++k) {
    //z = pcoord->x3v(k) - zc;
    for (int j=js; j<=je; ++j) {
      //y = pcoord->x2v(j) - yc;
      for (int i=is; i<=ie; ++i) {
        //x = pcoord->x1v(i) - xc;

        //r = std::sqrt( x*x + y*y + z*z );
        r = pcoord->x1v(i);

        if ( r < r0 ) {
          vr = r/t0;

          index = -1;
          mindv = 10.0e9;
          for (int l=0; l<NumToRead; l++) {
            dv = std::abs(vr-vr_in[l]);
            if ( dv < mindv ) {
              mindv = dv;
              index = l;
            }
          }
          printf("for r vr %5.3e %5.3e found neighbor id %d with rho vr %5.3e %5.3e\n",r,vr,index,rho_in[index],vr_in[index]);

          rho  = scaleFactorRho  * rho_in[index];
          temp = scaleFactorTemp * temp_in[index];

          if(!useGammaLaw) {
            mesaeos_dtget( &rho, &temp, &X, &Z, &use_solar, &fc12,
              &fn14, &fo16, &fne20, &presJunk, &Especific, &gamma);
            phydro->u(IEN,k,j,i) = Especific*rho + 0.5*rho*vr*vr;
          } else {
            // internal energy = ~10% kinetic
            phydro->u(IEN,k,j,i) = 1.1*0.5*rho*3.0e9*3.0e9;
          }

          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = rho*vr; // rho*x/t0;
          phydro->u(IM2,k,j,i) = 0.0; // rho*y/t0;
          phydro->u(IM3,k,j,i) = 0.0; // rho*z/t0;
        } else {
          phydro->u(IDN,k,j,i) = rhoFloor;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          if(!useGammaLaw) {
            phydro->u(IEN,k,j,i) = EspecificFloor*rhoFloor; // + 0.5*rho0R*vel0R*vel0R;
          } else {
            // internal energy < 10% kinetic
            phydro->u(IEN,k,j,i) = 1.1*0.5*rho*3.0e9*3.0e9;
          }
        }
      }
    }
  }
  return;
}
} // extern C
