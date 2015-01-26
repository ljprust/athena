//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "../mesh.hpp"

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // GetGamma
#include "../field/field.hpp"      // magnetic field

//======================================================================================
//! \file shock_tube.cpp
//  \brief Problem generator for shock tube problems.  
//
// Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
// shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  std::stringstream msg;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

// parse shock direction: {1,2,3} -> {x1,x2,x3}

  int shk_dir = pin->GetInteger("problem","shock_dir"); 

// parse shock location (must be inside grid)

  Real xshock = pin->GetReal("problem","xshock"); 
  if (shk_dir == 1 && (xshock < pmb->pmy_mesh->mesh_size.x1min ||
                       xshock > pmb->pmy_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 2 && (xshock < pmb->pmy_mesh->mesh_size.x2min ||
                       xshock > pmb->pmy_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (shk_dir == 3 && (xshock < pmb->pmy_mesh->mesh_size.x3min ||
                       xshock > pmb->pmy_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// Parse left state read from input file: dl,ul,vl,wl,[pl]

  Real wl[NFLUID+NFIELD];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) wl[IEN] = pin->GetReal("problem","pl");
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NFLUID  ] = pin->GetReal("problem","bxl");
    wl[NFLUID+1] = pin->GetReal("problem","byl");
    wl[NFLUID+2] = pin->GetReal("problem","bzl");
  }

// Parse right state read from input file: dr,ur,vr,wr,[pr]

  Real wr[NFLUID+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) wr[IEN] = pin->GetReal("problem","pr");
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NFLUID  ] = pin->GetReal("problem","bxr");
    wr[NFLUID+1] = pin->GetReal("problem","byr");
    wr[NFLUID+2] = pin->GetReal("problem","bzr");
  }

// Initialize the discontinuity in the Fluid variables ---------------------------------

  switch(shk_dir) {

//--- shock in 1-direction
  case 1:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pmb->x1v(i) < xshock) {
          pfl->u(IDN,k,j,i) = wl[IDN];
          pfl->u(IM1,k,j,i) = wl[IVX]*wl[IDN];
          pfl->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
          pfl->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wl[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        } else {
          pfl->u(IDN,k,j,i) = wr[IDN];
          pfl->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
          pfl->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
          pfl->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wr[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//--- shock in 2-direction
  case 2:
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (pmb->x2v(j) < xshock) {
        for (int i=is; i<=ie; ++i) {
          pfl->u(IDN,k,j,i) = wl[IDN];
          pfl->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
          pfl->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
          pfl->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wl[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }
      } else {
        for (int i=is; i<=ie; ++i) {
          pfl->u(IDN,k,j,i) = wr[IDN];
          pfl->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
          pfl->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
          pfl->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wr[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }
      }
    }}
    break;

//--- shock in 3-direction

  case 3:
    for (int k=ks; k<=ke; ++k) {
      if (pmb->x3v(k) < xshock) {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfl->u(IDN,k,j,i) = wl[IDN];
          pfl->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
          pfl->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
          pfl->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wl[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY] + wl[IVZ]*wl[IVZ]);
        }}
      } else {
        for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfl->u(IDN,k,j,i) = wr[IDN];
          pfl->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
          pfl->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
          pfl->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
          if (NON_BAROTROPIC_EOS) pfl->u(IEN,k,j,i) =
            wr[IEN]/(pfl->pf_eos->GetGamma() - 1.0)
            + 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY] + wr[IVZ]*wr[IVZ]);
        }}
      }
    }
    break;

  default:
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (shk_dir==1 && pmb->x1v(i) < xshock) {
          pfd->b.x1f(k,j,i) = wl[NFLUID  ];
          pfd->b.x2f(k,j,i) = wl[NFLUID+1];
          pfd->b.x3f(k,j,i) = wl[NFLUID+2];
        } else if (shk_dir==2 && pmb->x2v(j) < xshock) {
          pfd->b.x1f(k,j,i) = wl[NFLUID+2];
          pfd->b.x2f(k,j,i) = wl[NFLUID  ];
          pfd->b.x3f(k,j,i) = wl[NFLUID+1];
        } else if (shk_dir==3 && pmb->x3v(k) < xshock) {
          pfd->b.x1f(k,j,i) = wl[NFLUID+1];
          pfd->b.x2f(k,j,i) = wl[NFLUID+2];
          pfd->b.x3f(k,j,i) = wl[NFLUID];
        }

        if (shk_dir==1 && pmb->x1v(i) >= xshock) {
          pfd->b.x1f(k,j,i) = wr[NFLUID  ];
          pfd->b.x2f(k,j,i) = wr[NFLUID+1];
          pfd->b.x3f(k,j,i) = wr[NFLUID+2];
        } else if (shk_dir==2 && pmb->x2v(j) >= xshock) {
          pfd->b.x1f(k,j,i) = wr[NFLUID+2];
          pfd->b.x2f(k,j,i) = wr[NFLUID  ];
          pfd->b.x3f(k,j,i) = wr[NFLUID+1];
        } else if (shk_dir==3 && pmb->x3v(k) >= xshock)  {
          pfd->b.x1f(k,j,i) = wr[NFLUID+1];
          pfd->b.x2f(k,j,i) = wr[NFLUID+2];
          pfd->b.x3f(k,j,i) = wr[NFLUID];
        }
        if (NON_BAROTROPIC_EOS) {
          pfl->u(IEN,k,j,i) += 0.5*(SQR(pfd->b.x1f(k,j,i)) + SQR(pfd->b.x2f(k,j,i)) +
             SQR(pfd->b.x3f(k,j,i)));
        }
      }
    }}

// end by adding bi.x1 at ie+1, bi.x2 at je+1, and bi.x3 at ke+1

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pfd->b.x1f(k,j,ie+1) = pfd->b.x1f(k,j,ie);
    }}
    for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      pfd->b.x2f(k,je+1,i) = pfd->b.x2f(k,je,i);
    }}
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfd->b.x3f(ke+1,j,i) = pfd->b.x3f(ke,j,i);
    }}
  }

  return;
}
