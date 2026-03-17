//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file companion_gravity.cpp
//! \brief source terms due to gravity from a companion star

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../gravity/gravity.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::CompanionGravity
//! \brief Adds source terms for gravitational acceleration to conserved variables
//! \note
//! Assumes the companion is a point source on the polar axis.

void HydroSourceTerms::CompanionGravity(const Real dt,const AthenaArray<Real> *flux,
                                        const AthenaArray<Real> &prim,
                                        AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Gravity *pgrav = pmb->pgrav;

  // acceleration in x1-direction
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        Real r, theta, phi, z, x, y, xc, dist3, phi0;
        Real delPhi_x, delPhi_y, delPhi_z;
        Real delPhi_cyl, delPhi_r, delPhi_theta, delPhi_phi;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          r   = pmb->pcoord->x1v(i);
          phi = pmb->pcoord->x2v(j);
          z   = pmb->pcoord->x3v(k);
          x   = r*std::cos(phi);
          y   = r*std::sin(phi); 
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          r     = pmb->pcoord->x1v(i);
          theta = pmb->pcoord->x2v(j);
          phi   = pmb->pcoord->x3v(k);
          x     = r*std::sin(theta)*std::cos(phi);
          y     = r*std::sin(theta)*std::sin(phi);
          z     = r*std::cos(theta);
        }

        xc = r_companion;
        dist3 = std::pow((x-xc)*(x-xc)+y*y+z*z+r_plummer*r_plummer,1.5); 

        delPhi_x = gm_companion*(x-xc)/dist3;
        delPhi_y = gm_companion*y     /dist3; 
        delPhi_z = gm_companion*z     /dist3; 

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phi0 = delPhi_x*std::atan(delPhi_y/delPhi_x);
          delPhi_cyl = std::sqrt(delPhi_x*delPhi_x+delPhi_y*delPhi_y);
          delPhi_r = delPhi_cyl*std::cos(phi-phi0);          
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          delPhi_r = delPhi_x*std::sin(theta)*std::cos(phi)
                   + delPhi_y*std::sin(theta)*std::sin(phi)
                   + delPhi_z*std::cos(theta);
        }

        cons(IM1,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_r;
        cons(IEN,k,j,i) -= dt*0.5*(flux[X1DIR](IDN,k,j,i  )*delPhi_r + 
                                   flux[X1DIR](IDN,k,j,i+1)*delPhi_r);
      }
    }
  }

  if (pmb->block_size.nx2 > 1) {
    // acceleration in x2-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {

          Real r, theta, phi, z, x, y, xc, dist3, phi0;
          Real delPhi_x, delPhi_y, delPhi_z;
          Real delPhi_cyl, delPhi_r, delPhi_theta, delPhi_phi;

          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            r   = pmb->pcoord->x1v(i);
            phi = pmb->pcoord->x2v(j);
            z   = pmb->pcoord->x3v(k);
            x   = r*std::cos(phi);
            y   = r*std::sin(phi); 
          } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            r     = pmb->pcoord->x1v(i);
            theta = pmb->pcoord->x2v(j);
            phi   = pmb->pcoord->x3v(k);
            x     = r*std::sin(theta)*std::cos(phi);
            y     = r*std::sin(theta)*std::sin(phi);
            z     = r*std::cos(theta);
          }

          xc = r_companion;
          dist3 = std::pow((x-xc)*(x-xc)+y*y+z*z+r_plummer*r_plummer,1.5); 

          delPhi_x = gm_companion*(x-xc)/dist3;
          delPhi_y = gm_companion*y     /dist3; 
          delPhi_z = gm_companion*z     /dist3; 

          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            phi0 = delPhi_x*std::atan(delPhi_y/delPhi_x);
            delPhi_cyl = std::sqrt(delPhi_x*delPhi_x+delPhi_y*delPhi_y);
            delPhi_r = delPhi_cyl*std::cos(phi-phi0); 
            delPhi_phi = -delPhi_cyl*std::sin(phi-phi0);
            cons(IM2,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_phi;
            cons(IEN,k,j,i) -= dt*0.5*(flux[X2DIR](IDN,k,j  ,i)*delPhi_phi +
                                       flux[X2DIR](IDN,k,j+1,i)*delPhi_phi);
          } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            delPhi_theta = delPhi_x*std::cos(theta)*std::cos(phi)
                         + delPhi_y*std::cos(theta)*std::sin(phi)
                         - delPhi_z*std::sin(theta);
            cons(IM2,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_theta;
            cons(IEN,k,j,i) -= dt*0.5*(flux[X2DIR](IDN,k,j  ,i)*delPhi_theta +
                                       flux[X2DIR](IDN,k,j+1,i)*delPhi_theta);
          }
        }
      }
    }
  }

  if (pmb->block_size.nx3 > 1) {
    // acceleration in x3-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {

          Real r, theta, phi, z, x, y, xc, dist3, phi0;
          Real delPhi_x, delPhi_y, delPhi_z;
          Real delPhi_cyl, delPhi_r, delPhi_theta, delPhi_phi;

          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            r   = pmb->pcoord->x1v(i);
            phi = pmb->pcoord->x2v(j);
            z   = pmb->pcoord->x3v(k);
            x   = r*std::cos(phi);
            y   = r*std::sin(phi); 
          } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            r     = pmb->pcoord->x1v(i);
            theta = pmb->pcoord->x2v(j);
            phi   = pmb->pcoord->x3v(k);
            x     = r*std::sin(theta)*std::cos(phi);
            y     = r*std::sin(theta)*std::sin(phi);
            z     = r*std::cos(theta);
          }

          xc = r_companion;
          dist3 = std::pow((x-xc)*(x-xc)+y*y+z*z+r_plummer*r_plummer,1.5); 

          delPhi_x = gm_companion*(x-xc)/dist3;
          delPhi_y = gm_companion*y     /dist3; 
          delPhi_z = gm_companion*z     /dist3; 

          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            cons(IM3,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_z;
            cons(IEN,k,j,i) -= dt*0.5*(flux[X3DIR](IDN,k,j  ,i)*delPhi_z +
                                       flux[X3DIR](IDN,k+1,j,i)*delPhi_z);
          } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            delPhi_phi = -delPhi_x*std::sin(phi) + delPhi_y*std::cos(phi);
            cons(IM3,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_phi;
            cons(IEN,k,j,i) -= dt*0.5*(flux[X3DIR](IDN,k  ,j,i)*delPhi_phi +
                                       flux[X3DIR](IDN,k+1,j,i)*delPhi_phi);
          }
        }
      }
    }
  }

  return;
}
