//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shock_tube_bc.cpp
//! \brief Problem generator for shock tube problems.
//!
//! Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
//! shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt
#include <cstdio>
#include <iostream>   // endl
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
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// inflow/outflow BCs
void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
             Real time, Real dt,
             int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  Real rhoIn, vxIn, vyIn, vzIn, presIn;
} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  rhoIn  = pin->GetReal("problem","dl");
  vxIn   = pin->GetReal("problem","ul");
  vyIn   = pin->GetReal("problem","vl");
  vzIn   = pin->GetReal("problem","wl");
  presIn = pin->GetReal("problem","pl");

  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OuterX1);
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
        user_out_var(0,k,j,i) = pcoord->forceOnBoundary(IVX,k,j,i);
        user_out_var(1,k,j,i) = pcoord->forceOnBoundary(IVY,k,j,i);
        user_out_var(2,k,j,i) = pcoord->forceOnBoundary(IVZ,k,j,i);

        pcoord->forceOnBoundary(IDN,k,j,i) = 0.0;
        pcoord->forceOnBoundary(IVX,k,j,i) = 0.0;
        pcoord->forceOnBoundary(IVY,k,j,i) = 0.0;
        pcoord->forceOnBoundary(IVZ,k,j,i) = 0.0;
        pcoord->forceOnBoundary(IEN,k,j,i) = 0.0;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void OuterX1()
//  \brief Sets boundary condition on upstream boundary (oib) for wind tunnel
//
// Quantities at this boundary are held fixed at the constant upstream state

void OuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,iu+i) = rhoIn;
        prim(IM1,k,j,iu+i) = vxIn;
        prim(IM2,k,j,iu+i) = vyIn;
        prim(IM3,k,j,iu+i) = vzIn;
        prim(IEN,k,j,iu+i) = presIn;

      }
    }
  }
}

//========================================================================================
//! \fn Real press(Real rho, Real T)
//! \brief Calculate pressure as a function of density and temperature for H EOS.
//========================================================================================

Real press(Real rho, Real T) {
  // Ionization fraction
  Real x = 2. /(1 + std::sqrt(1 + 4. * rho * std::exp(1. / T) * std::pow(T, -1.5)));
  return rho * T * (1. + x);
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");

  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");

  // Parse left state read from input file: dl,ul,vl,wl,[pl]
  Real wl[NHYDRO+NFIELD];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) {
    if (pin->DoesParameterExist("problem","Tl"))
      wl[IPR] = press(wl[IDN], pin->GetReal("problem","Tl"));
    else
      wl[IPR] = pin->GetReal("problem","pl");
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NHYDRO  ] = pin->GetReal("problem","bxl");
    wl[NHYDRO+1] = pin->GetReal("problem","byl");
    wl[NHYDRO+2] = pin->GetReal("problem","bzl");
  }

  // Parse right state read from input file: dr,ur,vr,wr,[pr]
  Real wr[NHYDRO+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) {
    if (pin->DoesParameterExist("problem","Tr"))
      wr[IPR] = press(wr[IDN], pin->GetReal("problem","Tr"));
    else
      wr[IPR] = pin->GetReal("problem","pr");
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NHYDRO  ] = pin->GetReal("problem","bxr");
    wr[NHYDRO+1] = pin->GetReal("problem","byr");
    wr[NHYDRO+2] = pin->GetReal("problem","bzr");
  }

  // Initialize the discontinuity in the Hydro variables ---------------------------------

  switch(shk_dir) {
    //--- shock in 1-direction
    case 1:
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (pcoord->x1v(i) < xshock) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            } else {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;
      //--- shock in 2-direction
    case 2:
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          if (pcoord->x2v(j) < xshock) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            }
          } else {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;

      //--- shock in 3-direction
    case 3:
      for (int k=ks; k<=ke; ++k) {
        if (pcoord->x3v(k) < xshock) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            }
          }
        } else {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;

    default:
      msg << "### FATAL ERROR in Problem Generator" << std::endl
          << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
      ATHENA_ERROR(msg);
  }

  // uniformly fill all scalars to have equal concentration
  // mass fraction? or concentration?
  constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;
  if (NSCALARS > 0) {
    for (int n=0; n<NSCALARS; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pscalars->s(n,k,j,i) = 1.0/scalar_norm*phydro->u(IDN,k,j,i);
          }
        }
      }
    }
  }
  return;
}
