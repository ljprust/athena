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
#include <algorithm>
#include <cmath>
#include <cstring>    // memset
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/fft_gravity.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// inflow/outflow BCs
void WindTunnel2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void WindTunnel2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void WindTunnel2DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void WindTunnel2DInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void WindTunnel2DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void WindTunnel2DInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                       Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void PointExplode(MeshBlock *pmb, const Real time, const Real dt,
                       const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                       const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                       AthenaArray<Real> &cons_scalar);
void VacuumSource(MeshBlock *pmb, const Real time, const Real dt,
                       const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                       const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                       AthenaArray<Real> &cons_scalar);

void StaticInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
// problem parameters which are useful to make global to this file
Real gm0, rho0, vel0, p0, gammagas, semimajor, gmstar;
bool diode, hydrostatic, staticBoundary, expgrad;
Real pvacuum, dvacuum, densgrad;
Real xvac, yvac, zvac, rvac, gm_point;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real four_pi_G = pin->GetReal("problem","four_pi_G");
  SetFourPiG(four_pi_G);
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
  //EnrollUserBoundaryFunction(BoundaryFace::outer_x1, WindTunnel2DOuterX1);
  //EnrollUserBoundaryFunction(BoundaryFace::inner_x1, WindTunnel2DInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, StaticInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, StaticOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, StaticInnerX2);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x2, StaticOuterX2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, StaticInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, StaticOuterX3);
  xvac = pin->GetOrAddReal("problem","xvac",0.0);
  yvac = pin->GetOrAddReal("problem","yvac",0.0);
  zvac = pin->GetOrAddReal("problem","zvac",0.0);
  rvac = pin->GetOrAddReal("problem","rvac",0.0);
  gm_point = pin->GetOrAddReal("problem","gm_point",0.0);
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
//  \brief Sets boundary condition on upstream boundary for wind tunnel
//
// Quantities at this boundary are held fixed at the constant upstream state

void WindTunnel2DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

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

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DInnerX1()
//  \brief Sets boundary condition on upstream boundary for wind tunnel
//
// Quantities at this boundary are held fixed at the constant upstream state

void WindTunnel2DInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rho, pres;
  Real cs2, Minf2, ratio, y;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          y = pco->x2v(j);
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "invalid coordinate system" << std::endl;
          ATHENA_ERROR(msg);
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
          prim(IDN,k,j,il-i) = rho;
          prim(IM1,k,j,il-i) = vel0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = 0.0;
          prim(IEN,k,j,il-i) = pres;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DOuterX2()
//  \brief Sets boundary condition on outflow boundary for wind tunnel

void WindTunnel2DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {

        prim(IDN,k,ju+j,i) = prim(IDN,k,ju,i);
        prim(IM1,k,ju+j,i) = prim(IM1,k,ju,i);
        prim(IM3,k,ju+j,i) = prim(IM3,k,ju,i);
        prim(IEN,k,ju+j,i) = prim(IEN,k,ju,i);

        // ensure that no gas enters through the outflow boundary
        // by giving the radial velocity some TLC
        applyDiode = prim(IM2,k,ju,i) < 0.0;
        if (diode && applyDiode) {
          prim(IM2,k,ju+j,i) = 0.0;
        } else {
          prim(IM2,k,ju+j,i) = prim(IM2,k,ju,i);
        }

      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DInnerX2()
//  \brief Sets boundary condition on outflow boundary for wind tunnel

void WindTunnel2DInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {

        prim(IDN,k,jl-j,i) = prim(IDN,k,jl,i);
        prim(IM1,k,jl-j,i) = prim(IM1,k,jl,i);
        prim(IM3,k,jl-j,i) = prim(IM3,k,jl,i);
        prim(IEN,k,jl-j,i) = prim(IEN,k,jl,i);

        // ensure that no gas enters through the outflow boundary
        // by giving the radial velocity some TLC
        applyDiode = prim(IM2,k,jl,i) > 0.0;
        if (diode && applyDiode) {
          prim(IM2,k,jl-j,i) = 0.0;
        } else {
          prim(IM2,k,jl-j,i) = prim(IM2,k,jl,i);
        }

      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DOuterX3()
//  \brief Sets boundary condition on outflow boundary for wind tunnel

void WindTunnel2DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {

        prim(IDN,ku+k,j,i) = prim(IDN,ku,j,i);
        prim(IM1,ku+k,j,i) = prim(IM1,ku,j,i);
        prim(IM2,ku+k,j,i) = prim(IM2,ku,j,i);
        prim(IEN,ku+k,j,i) = prim(IEN,ku,j,i);

        // ensure that no gas enters through the outflow boundary
        // by giving the radial velocity some TLC
        applyDiode = prim(IM3,ku,j,i) < 0.0;
        if (diode && applyDiode) {
          prim(IM3,ku+k,j,i) = 0.0;
        } else {
          prim(IM3,ku+k,j,i) = prim(IM3,ku,j,i);
        }

      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void WindTunnel2DInnerX3()
//  \brief Sets boundary condition on outflow boundary for wind tunnel

void WindTunnel2DInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {

        prim(IDN,kl-k,j,i) = prim(IDN,kl,j,i);
        prim(IM1,kl-k,j,i) = prim(IM1,kl,j,i);
        prim(IM2,kl-k,j,i) = prim(IM2,kl,j,i);
        prim(IEN,kl-k,j,i) = prim(IEN,kl,j,i);

        // ensure that no gas enters through the outflow boundary
        // by giving the radial velocity some TLC
        applyDiode = prim(IM3,kl,j,i) > 0.0;
        if (diode && applyDiode) {
          prim(IM3,kl-k,j,i) = 0.0;
        } else {
          prim(IM3,kl-k,j,i) = prim(IM3,kl,j,i);
        }

      }
    }
  }
}

void PointExplode(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

  Real rad, theta, phi;
  Real x, y, z;
  Real r2_relative;
  bool inExp;
  Real tend;
  Real rho, v1, v2, v3, vsq;

  tend = tstartexp + dtexp;
  bool ininterval = time > tstartexp && time < tend;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          rad = pmb->pcoord->x1v(i);
          theta = pmb->pcoord->x2v(j);
          phi = pmb->pcoord->x3v(k);
          x = rad*std::sin(theta)*std::cos(phi);
          y = rad*std::sin(theta)*std::sin(phi);
          z = rad*std::cos(theta);
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "invalid coordinate system" << std::endl;
          ATHENA_ERROR(msg);
        }
        x = x-dexp;
        r2_relative = x*x + y*y + z*z;
        inExp = rexp*rexp > r2_relative;
        rho = prim(IDN,k,j,i);
        v1  = prim(IVX,k,j,i);
        v2  = prim(IVY,k,j,i);
        v3  = prim(IVZ,k,j,i);
        vsq = v1*v1 + v2*v2 + v3*v3;
        if (NON_BAROTROPIC_EOS && Pexp>0.0 && inExp && ininterval) {
          //cons(IEN,k,j,i) -= dt*den*SQR(Omega_0)*prim(IVZ,k,j,i)*x3*fsmooth;
          cons(IEN,k,j,i) = Pexp/(gammagas-1.0) + 0.5*rho*vsq;
        }
      }
    }
  }
  return;
}

void VacuumSource(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

  Real x, y, z;
  Real r2_relative;
  bool invac;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pmb->pcoord->x1v(i) - xvac;
          y = pmb->pcoord->x2v(j) - yvac;
          z = pmb->pcoord->x3v(k) - zvac;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in windtunnel.cpp ProblemGenerator" << std::endl
              << "invalid coordinate system" << std::endl;
          ATHENA_ERROR(msg);
        }
        r2_relative = x*x + y*y + z*z;
        invac = rvac*rvac > r2_relative;
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

void StaticInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        prim(IDN,k,j,il-i) = rho0;
        prim(IM1,k,j,il-i) = 0.0;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = p0;
      }
    }
  }
}

void StaticOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        prim(IDN,k,j,iu+i) = rho0;
        prim(IM1,k,j,iu+i) = 0.0;
        prim(IM2,k,j,iu+i) = 0.0;
        prim(IM3,k,j,iu+i) = 0.0;
        prim(IEN,k,j,iu+i) = p0;
      }
    }
  }
}

void StaticInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,jl-j,i) = rho0;
        prim(IM1,k,jl-j,i) = 0.0;
        prim(IM2,k,jl-j,i) = 0.0;
        prim(IM3,k,jl-j,i) = 0.0;
        prim(IEN,k,jl-j,i) = p0;
      }
    }
  }
}

void StaticOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,ju+j,i) = rho0;
        prim(IM1,k,ju+j,i) = 0.0;
        prim(IM2,k,ju+j,i) = 0.0;
        prim(IM3,k,ju+j,i) = 0.0;
        prim(IEN,k,ju+j,i) = p0;
      }
    }
  }
}

void StaticInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,kl-k,j,i) = rho0;
        prim(IM1,kl-k,j,i) = 0.0;
        prim(IM2,kl-k,j,i) = 0.0;
        prim(IM3,kl-k,j,i) = 0.0;
        prim(IEN,kl-k,j,i) = p0;
      }
    }
  }
}

void StaticOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,ku+k,j,i) = rho0;
        prim(IM1,ku+k,j,i) = 0.0;
        prim(IM2,ku+k,j,i) = 0.0;
        prim(IM3,ku+k,j,i) = 0.0;
        prim(IEN,ku+k,j,i) = p0;
      }
    }
  } 
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
