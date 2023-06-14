//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file user_gravity.cpp
//! \brief source terms due to user-defined gravity

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
//! \fn void HydroSourceTerms::UserGravity
//! \brief Adds source terms for user-defined gravitational acceleration to conserved variables
//! \note
//! Currently set up for gravity from a point source in the far-field
//! where the center of the domain is a distance /semimajor/ from the point
//! source, which has GM = gmstar. Assumes the direction of gravity is constant
//! but not the magnitude.

void HydroSourceTerms::UserGravity(const Real dt,const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Gravity *pgrav = pmb->pgrav;
  Real r, x, y, z, delPhix, delPhiy, delPhiz;

  // acceleration in 1-direction
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        x = pmb->pcoord->x1v(i) - xvac;
        y = pmb->pcoord->x2v(j) - yvac;
        z = pmb->pcoord->x3v(k) - zvac;
        r = std::sqrt(x*x + y*y + z*z);
        delPhix = -gm_point/r/r/r*x;
        delPhiy = -gm_point/r/r/r*y;
        delPhiz = -gm_point/r/r/r*z;

        cons(IM1,k,j,i) += dt*prim(IDN,k,j,i)*delPhix;
        cons(IM2,k,j,i) += dt*prim(IDN,k,j,i)*delPhiy;
        cons(IM3,k,j,i) += dt*prim(IDN,k,j,i)*delPhiz;

        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += dt*0.5*(delPhix*(flux[X1DIR](IDN,k,j,i  )  +
                                              flux[X1DIR](IDN,k,j,i+1)) +
                                     delPhiy*(flux[X2DIR](IDN,k,j,i  )  +
                                              flux[X2DIR](IDN,k,j+1,i)) +
                                     delPhiz*(flux[X3DIR](IDN,k,j,i  )  +
                                              flux[X3DIR](IDN,k+1,j,i)));

      }
    }
  }

  return;
}
/*
void HydroSourceTerms::UserGravity(const Real dt,const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Gravity *pgrav = pmb->pgrav;

  // acceleration in 1-direction
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real dx1 = pmb->pcoord->dx1v(i);
        Real dtodx1 = dt/dx1;

        // define phi for center and neighbors
        Real r, theta, phi, y;
        r = pmb->pcoord->x1v(i);
        theta = pmb->pcoord->x2v(j);
        phi = pmb->pcoord->x3v(k);
        y = r*std::sin(theta)*std::cos(phi);

        Real delPhi = gmstar/SQR(semimajor+y);
        Real delPhi_r = delPhi*std::sin(theta)*std::cos(phi);
        cons(IM1,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_r;
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) -= dt*0.5*(flux[X1DIR](IDN,k,j,i  )*delPhi_r +
                                     flux[X1DIR](IDN,k,j,i+1)*delPhi_r);

        // Real phic = -gmstar / (semimajor+y);
        //
        // Real rL, yL;
        // rL = pmb->pcoord->x1v(i-1);
        // yL = rL*std::sin(theta)*std::cos(phi);
        // Real phiLeft = -gmstar / (semimajor+yL);
        //
        // Real rR, yR;
        // rR = pmb->pcoord->x1v(i+1);
        // yR = rR*std::sin(theta)*std::cos(phi);
        // Real phiRight = -gmstar / (semimajor+yR);
        //
        // Real phil = 0.5*(phiLeft+phic);
        // Real phir = 0.5*(phic+phiRight);
        // cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phir-phil);
        // //std::cout << dtodx1*prim(IDN,k,j,i)*(phir-phil) << std::endl;
        //
        // if (NON_BAROTROPIC_EOS)
        //   cons(IEN,k,j,i) -= dtodx1*(flux[X1DIR](IDN,k,j,i  )*(phic - phil) +
        //                              flux[X1DIR](IDN,k,j,i+1)*(phir - phic));

      }
    }
  }

  if (pmb->block_size.nx2 > 1) {
    // acceleration in 2-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx2 = pmb->pcoord->dx2v(j);
          Real dtodx2 = dt/dx2;

          // define phi for center and neighbors
          Real r, theta, phi, y;
          r = pmb->pcoord->x1v(i);
          theta = pmb->pcoord->x2v(j);
          phi = pmb->pcoord->x3v(k);
          y = r*std::sin(theta)*std::cos(phi);
          //Real phic = -gmstar / (semimajor+y);

          Real delPhi = gmstar/SQR(semimajor+y);
          Real delPhi_theta = delPhi*std::cos(theta)*std::cos(phi);
          cons(IM2,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_theta;
          if (NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) -= dt*0.5*(flux[X2DIR](IDN,k,j  ,i)*delPhi_theta +
                                       flux[X2DIR](IDN,k,j+1,i)*delPhi_theta);

          // Real thetaL, yL;
          // thetaL = pmb->pcoord->x2v(j-1);
          // yL = r*std::sin(thetaL)*std::cos(phi);
          // Real phiLeft = -gmstar / (semimajor+yL);
          //
          // Real thetaR, yR;
          // thetaR = pmb->pcoord->x2v(j+1);
          // yR = r*std::sin(thetaR)*std::cos(phi);
          // Real phiRight = -gmstar / (semimajor+yR);
          //
          // Real phil = 0.5*(phiLeft+phic);
          // Real phir = 0.5*(phic+phiRight);
          // cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phir-phil);
          // if (NON_BAROTROPIC_EOS)
          //   cons(IEN,k,j,i) -= dtodx2*(flux[X2DIR](IDN,k,j  ,i)*(phic - phil) +
          //                              flux[X2DIR](IDN,k,j+1,i)*(phir - phic));

        }
      }
    }
  }

  if (pmb->block_size.nx3 > 1) {
    // acceleration in 3-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx3 = pmb->pcoord->dx3v(k);
          Real dtodx3 = dt/dx3;

          // define phi for center and neighbors
          Real r, theta, phi, y;
          r = pmb->pcoord->x1v(i);
          theta = pmb->pcoord->x2v(j);
          phi = pmb->pcoord->x3v(k);
          y = r*std::sin(theta)*std::cos(phi);
          //Real phic = -gmstar / (semimajor+y);

          Real delPhi = gmstar/SQR(semimajor+y);
          Real delPhi_phi = -delPhi*std::sin(phi);
          cons(IM3,k,j,i) -= dt*prim(IDN,k,j,i)*delPhi_phi;
          if (NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) -= dt*0.5*(flux[X3DIR](IDN,k  ,j,i)*delPhi_phi +
                                       flux[X3DIR](IDN,k+1,j,i)*delPhi_phi);


          // Real phiL, yL;
          // phiL = pmb->pcoord->x3v(k-1);
          // yL = r*std::sin(theta)*std::cos(phiL);
          // Real phiLeft = -gmstar / (semimajor+yL);
          //
          // Real phiR, yR;
          // phiR = pmb->pcoord->x3v(k+1);
          // yR = r*std::sin(theta)*std::cos(phiR);
          // Real phiRight = -gmstar / (semimajor+yR);
          //
          // Real phil = 0.5*(phiLeft+phic);
          // Real phir = 0.5*(phic+phiRight);
          // cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phir-phil);
          // if (NON_BAROTROPIC_EOS)
          //   cons(IEN,k,j,i) -= dtodx3*(flux[X3DIR](IDN,k  ,j,i)*(phic - phil) +
          //                              flux[X3DIR](IDN,k+1,j,i)*(phir - phic));

        }
      }
    }
  }

  return;
}
*/
