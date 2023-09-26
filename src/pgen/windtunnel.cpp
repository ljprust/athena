//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file windtunnel.cpp
//! \brief Initializes parallel flow in one direction in both cylindrical and
//! spherical polar coordinates.

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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// inflow/outflow BCs
void WindTunnel2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void WindTunnel2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void VacuumSource(MeshBlock *pmb, const Real time, const Real dt,
                       const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                       const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                       AthenaArray<Real> &cons_scalar);

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
// problem parameters which are useful to make global to this file
Real gm0, rho0, vel0, p0, gammagas, semimajor, gmstar;
bool diode, hydrostatic, staticBoundary, expgrad;
Real pvacuum, dvacuum, densgrad, rvacuum;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  rho0 = pin->GetOrAddReal("problem","rho0",1.0);
  vel0 = pin->GetOrAddReal("problem","vel0",1.0);
  p0 = pin->GetOrAddReal("problem","p0",1.0);
  gammagas = pin->GetOrAddReal("hydro","gamma",0.0);
  diode = pin->GetOrAddBoolean("problem","diode",false);
  hydrostatic = pin->GetOrAddBoolean("problem","hydrostatic",false);
  semimajor = pin->GetOrAddReal("problem","semimajor",0.0);
  gmstar = pin->GetOrAddReal("problem","gm_star",0.0);
  pvacuum = pin->GetOrAddReal("problem","pvacuum",0.0);
  dvacuum = pin->GetOrAddReal("problem","dvacuum",0.0);
  densgrad = pin->GetOrAddReal("problem","densgrad",0.0);
  expgrad = pin->GetOrAddBoolean("problem","expgrad",false);
  staticBoundary = pin->GetOrAddBoolean("problem","staticBoundary",false);
  rvacuum = pin->GetOrAddReal("problem","rvacuum",0.0);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, WindTunnel2DOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, WindTunnel2DInnerX1);
  EnrollUserExplicitSourceFunction(VacuumSource);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes wind tunnel.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real x1, x2, x3;
  Real rho, pres;
  Real cs2, Minf2, ratio, y;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          y = x1*std::sin(x2);
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          y = x1*std::sin(x2)*std::cos(x3);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          y = x2;
        }

        if (hydrostatic) {
          cs2 = gammagas*p0/rho0;
          Minf2 = gmstar/semimajor/cs2; // vel0*vel0/cs2;
          ratio = 1.0 - (gammagas-1.0)*Minf2/(1.0+semimajor/y);
          rho = rho0*std::pow( ratio, 1.0/(gammagas-1.0) );
          pres = p0*std::pow( ratio, gammagas/(gammagas-1.0) );
        } else if (expgrad) {
          rho = rho0*std::exp(densgrad*y);
          pres = p0*std::pow( rho0, gammagas-1.0 )*std::exp(densgrad*y*gammagas);
        } else {
          rho = rho0*(1.0-y*densgrad);
          pres = p0;
        }

        phydro->u(IDN,k,j,i) = rho;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          //GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
          phydro->u(IM1,k,j,i) =  rho*vel0*std::cos(x2); // radial
          phydro->u(IM2,k,j,i) = -rho*vel0*std::sin(x2); // azimuth
          phydro->u(IM3,k,j,i) =  0.0;               // z
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM1,k,j,i) =  rho*vel0*std::cos(x2); // radial
          phydro->u(IM2,k,j,i) = -rho*vel0*std::sin(x2); // polar
          phydro->u(IM3,k,j,i) =  0.0;               // azimuth
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          phydro->u(IM1,k,j,i) =  rho*vel0; // x
          phydro->u(IM2,k,j,i) =  0.0; // y
          phydro->u(IM3,k,j,i) =  0.0; // z
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "no acceptable coord system found" << std::endl;
          ATHENA_ERROR(msg);
        }

        phydro->u(IEN,k,j,i) = pres/(gammagas-1.0) + 0.5*rho*vel0*vel0;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DOuterX1()
//  \brief Sets boundary condition on upstream boundary (oib) for wind tunnel
//
// Quantities at this boundary are held fixed at the constant upstream state

void WindTunnel2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool inflow, applyDiode;
  Real phi;
  Real rho, pres;
  Real cs2, Minf2, ratio, y;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        //rad=pco->x1v(iu+i);
        phi=pco->x2v(j);
        //z=pco->x3v(k);
        inflow = phi > 3.14159/2.0 && phi < 3.14159*3.0/2.0;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          y = pco->x2v(iu+i)*std::sin(pco->x1v(j));
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          y = pco->x1v(iu+i)*std::sin(pco->x2v(j))*std::cos(pco->x3v(k));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          y = pco->x2v(j);
        }

        if (hydrostatic) {
          cs2 = gammagas*p0/rho0;
          Minf2 = gmstar/semimajor/cs2; // vel0*vel0/cs2;
          ratio = 1.0 - (gammagas-1.0)*Minf2/(1.0+semimajor/y);
          rho = rho0*std::pow( ratio, 1.0/(gammagas-1.0) );
          pres = p0*std::pow( ratio, gammagas/(gammagas-1.0) );
        } else if (expgrad) {
          rho = rho0*std::exp(densgrad*y);
          pres = p0*std::pow( rho0, gammagas-1.0 )*std::exp(densgrad*y*gammagas);
        } else {
          rho = rho0*(1.0-y*densgrad);
          pres = p0;
        }

        if (inflow || staticBoundary) {
          // set half of the outer boundary to upstream conditions
          prim(IDN,k,j,iu+i) =  rho;
          prim(IM1,k,j,iu+i) =  vel0*std::cos(phi); // radial
          prim(IM2,k,j,iu+i) = -vel0*std::sin(phi); // azimuth
          prim(IM3,k,j,iu+i) =  0.0;               // z
          prim(IEN,k,j,iu+i) =  pres;
        } else {
          // the other half lets gas outflow freely
          prim(IDN,k,j,iu+i) = prim(IDN,k,j,iu);
          prim(IM2,k,j,iu+i) = prim(IM2,k,j,iu);
          prim(IM3,k,j,iu+i) = prim(IM3,k,j,iu);
          prim(IEN,k,j,iu+i) = prim(IEN,k,j,iu);

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
//! \fn void WindTunnel2DInnerX1()
//  \brief Sets vacuum inner boundary
//
// Quantities in ghost cells are set to some pressure and density

void WindTunnel2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rho, pres;

  if ( pvacuum==0.0 || dvacuum==0.0 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
        << "vacuum pressure and/or density not set" << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = dvacuum;
        prim(IM1,k,j,il-i) = 0.0;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pvacuum;

      }
    }
  }
}

void VacuumSource(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

  Real r;
  bool invac;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          r = pmb->pcoord->x1v(i);
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "invalid coordinate system" << std::endl;
          ATHENA_ERROR(msg);
        }
        invac = rvacuum > r;
        if (invac) {
          cons(IDN,k,j,i) = dvacuum;
          cons(IVX,k,j,i) = 0.0;
          cons(IVY,k,j,i) = 0.0;
          cons(IVZ,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) = pvacuum/(gammagas-1.0);
          }
        }
      }
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

}
