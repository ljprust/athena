//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ablation.cpp
//! \brief Initializes parallel flow in one direction in spherical polar coordinates 
//         with planetary material on surface of hard boundary.

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

// inflow/outflow BCs
void WindTunnelInflowOutflow(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                             Real time, Real dt,
                             int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void PlanetInnerBoundary(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                         Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// store planet profile in these
std::vector<Real> r_in, rho_in, pres_in;

namespace {
// problem parameters which are useful to make global to this file
Real rho_inf, vel_inf, pres_inf, gammagas; 
Real R_planet, GM_planet, rho_inner, pres_inner;
bool diode, staticBoundary;
int  NumToRead;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // read in parameters from parameter file
  GM_planet      = pin->GetOrAddReal(   "problem", "GM",             0.0);
  rho_inf        = pin->GetOrAddReal(   "problem", "rho_inf",        1.0);
  vel_inf        = pin->GetOrAddReal(   "problem", "vel_inf",        1.0);
  pres_inf       = pin->GetOrAddReal(   "problem", "pres_inf",       1.0);
  gammagas       = pin->GetOrAddReal(   "hydro",   "gamma",          0.0);
  diode          = pin->GetOrAddBoolean("problem", "diode",          false);
  staticBoundary = pin->GetOrAddBoolean("problem", "staticBoundary", false);
  rho_inner      = pin->GetOrAddReal(   "problem", "rho_inner",      0.0);
  pres_inner     = pin->GetOrAddReal(   "problem", "pres_inner",     0.0);
  NumToRead      = pin->GetOrAddInteger("problem", "NumToRead",      0);

  // enroll boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, WindTunnelInflowOutflow);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, PlanetInnerBoundary);

  // read in planet profile from text files
  char rFile[256], rhoFile[256], presFile[256];
  sprintf(rFile,    "r.txt");
  sprintf(rhoFile,  "rho.txt");
  sprintf(presFile, "pres.txt");
  printf("Opening data files with state variables...\n");
  std::ifstream rFileRead, rhoFileRead, presFileRead;
  rFileRead.open(rFile);
  rhoFileRead.open(rhoFile);
  presFileRead.open(presFile);
  Real r, rho, pres;
  for (int l = 0; l < NumToRead; l++) {
    rFileRead    >> r;
    rhoFileRead  >> rho;
    presFileRead >> pres;
    r_in.push_back(r);
    rho_in.push_back(rho);
    pres_in.push_back(pres);
  }
  printf("Done reading, closing data files\n");
  rFileRead.close();
  rhoFileRead.close();
  presFileRead.close();

  // get the radius of the planet surface
  R_planet = r_in[NumToRead-1];

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes wind tunnel.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x1, x2, x3;
  Real rho, pres, vel, scalar;
  int index;
  Real r_L, r_R, rho_L, rho_R, pres_L, pres_R;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);

        if (x1 < R_planet) { // planetary material

          // find our location in the planet profile
          index = -1;
          for (int l=0; l<NumToRead; ++l) {
            if (r_in[l]>x1) {
              index = l;
              break;
            }
          }

          // interpolate the fluid variables
          r_L    = r_in[index-1];
          r_R    = r_in[index-1];
          rho_L  = rho_in[index-1];
          rho_R  = rho_in[index];
          pres_L = pres_in[index-1];
          pres_R = pres_in[index];
          rho    = rho_L  + (rho_R-rho_L)   * (x1-r_L) / (r_R-r_L);
          pres   = pres_L + (pres_R-pres_L) * (x1-r_L) / (r_R-r_L);

          vel    = 0.0; // static
          scalar = 1.0; // planetary material

        } else { // disk material

          rho    = rho_inf;
          pres   = pres_inf;
          vel    = vel_inf;
          scalar = 0.0; // not planetary material

        }

        // finally, set all fluid variables
        phydro->u(IDN,k,j,i) = rho;
        pscalars->s(0,k,j,i) = scalar*rho;
        phydro->u(IEN,k,j,i) = pres/(gammagas-1.0) + 0.5*rho*vel*vel;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM1,k,j,i) =  rho*vel*std::cos(x2); // radial
          phydro->u(IM2,k,j,i) = -rho*vel*std::sin(x2); // azimuth
          phydro->u(IM3,k,j,i) =  0.0;               // z
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM1,k,j,i) =  rho*vel*std::cos(x2); // radial
          phydro->u(IM2,k,j,i) = -rho*vel*std::sin(x2); // polar
          phydro->u(IM3,k,j,i) =  0.0;               // azimuth
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          phydro->u(IM1,k,j,i) =  rho*vel; // x
          phydro->u(IM2,k,j,i) =  0.0; // y
          phydro->u(IM3,k,j,i) =  0.0; // z
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in ablation.cpp ProblemGenerator" << std::endl
              << "no acceptable coord system found" << std::endl;
          ATHENA_ERROR(msg);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnelInflowOutflow()
//  \brief Sets upstream and downstream boundary conditions

void WindTunnelInflowOutflow(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                             Real time, Real dt,
                             int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool inflow, applyDiode;
  Real theta;
  Real rho, pres;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        //rad = pco->x1v(iu+i);
        theta = pco->x2v(j);

        // set half of boundary to inflow
        inflow = theta > 3.14159/2.0;

        if (inflow || staticBoundary) {
          // set half of the outer boundary to upstream conditions
          prim(IDN,k,j,iu+i) =  rho_inf;
          prim(IM1,k,j,iu+i) =  vel_inf*std::cos(theta);
          prim(IM2,k,j,iu+i) = -vel_inf*std::sin(theta);
          prim(IM3,k,j,iu+i) =  0.0;
          prim(IEN,k,j,iu+i) =  pres_inf;
          for (int l=0; l<NSCALARS; ++l) {
            pmb->pscalars->r(l,k,j,iu+i) = 0.0;
          }
        } else {
          // the other half lets gas outflow freely
          prim(IDN,k,j,iu+i) = prim(IDN,k,j,iu);
          prim(IM2,k,j,iu+i) = prim(IM2,k,j,iu);
          prim(IM3,k,j,iu+i) = prim(IM3,k,j,iu);
          prim(IEN,k,j,iu+i) = prim(IEN,k,j,iu);
          for (int l=0; l<NSCALARS; ++l) {
            pmb->pscalars->r(l,k,j,iu+i) = pmb->pscalars->r(l,k,j,iu);
          }

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
}

//----------------------------------------------------------------------------------------
//! \fn void PlanetInnerBoundary()
//  \brief Sets constant inner boundary

void PlanetInnerBoundary(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                         Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = rho_inner;
        prim(IM1,k,j,il-i) = 0.0;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pres_inner;

        pmb->pscalars->r(0,k,j,i) = 1.0;

      }
    }
  }
}

//namespace {
//}
