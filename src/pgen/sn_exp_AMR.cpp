//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sn.cpp
//! \brief Initializes SN ejecta into polar wedge

// C headers

// C++ headers
#include <algorithm>  // min
#include <cstdio>
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

// outflow condition with diode
void SNOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// inject SN ejecta
void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// set time step
Real MyTimeStep(MeshBlock* pmb);

//set refinement condition
int RefinementCondition(MeshBlock *pmb);

namespace {
// problem parameters which are useful to make global to this file
Real rho0, p0, gammagas, gasEntropy;
bool diode;
Real Rej, Eej, Mej, vmax;
Real smallDt, largeDt, smallDtDuration;

// new high-velocity tail fit parameters
Real v0_tail, vc_tail, t0_tail, rho0_tail;

// derived normalization A = rho0_tail * t0_tail^3
Real rho_norm_tail;

// minimum physical time (seconds) used in the fit
Real t_phys_floor;

// physical time at simulation t=0 (seconds)
Real t_phys0; 

constexpr int I_EJ = 0; //ejecta tracer

Real refine_drho_frac;
Real refine_dej_jump;
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

// read parameters from input file
void Mesh::InitUserMeshData(ParameterInput *pin) {

 
  // ambient medium properties
  rho0       = pin->GetOrAddReal("problem","rho0",1.0);
  p0         = pin->GetOrAddReal("problem","p0",1.0);

  // information for thermal energy
  gammagas   = pin->GetReal("hydro","gamma");
  gasEntropy = pin->GetReal("problem","gasEntropy");

  // use diode condition for outflow
  diode      = pin->GetOrAddBoolean("problem","diode",false);

  // properties of ejecta
  Rej        = pin->GetOrAddReal("problem","Rejecta",1.0);
  Eej        = pin->GetOrAddReal("problem","Eejecta",1.0);
  Mej        = pin->GetOrAddReal("problem","Mejecta",1.0);
  vmax       = pin->GetOrAddReal("problem","vmax",1.0);

  // NEW: tail-fit parameters (your v0, vc, t0, rho0)
  v0_tail    = pin->GetOrAddReal("problem","v0_tail",4.2e9);
  vc_tail    = pin->GetOrAddReal("problem","vc_tail",4.5e8);
  t0_tail    = pin->GetOrAddReal("problem","t0_tail",73.0);
  rho0_tail  = pin->GetOrAddReal("problem","rho0_tail",5.2e-5);

  // earliest physical time represented in the fit (seconds)
  t_phys_floor = pin->GetOrAddReal("problem","t_phys_floor",10.0);
  t_phys0 = pin->GetOrAddReal("problem","t_phys0", Rej / vmax);

  // prefactor A = rho0_tail * t0_tail^3, used as rho_norm_tail * t^{-3}
  rho_norm_tail = rho0_tail * std::pow(t0_tail, 3);

  // properties of initial timestep
  smallDt = pin->GetReal("problem", "smallDt");
  largeDt = pin->GetReal("problem", "largeDt");
  smallDtDuration = pin->GetReal("problem", "smallDtDuration");

  //refinement parameters
  refine_drho_frac = pin->GetOrAddReal("problem","refine_drho_frac",0.1);
  refine_dej_jump  = pin->GetOrAddReal("problem","refine_dej_jump",0.4);

  // tell athena to use our user-defined boundaries
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);

  // set a very small timestep initially
  EnrollUserTimeStepFunction(MyTimeStep);

  //enroll refinement condition
  EnrollUserRefinementCondition(RefinementCondition);


  if (Globals::my_rank == 0) {
    std::printf("### DEBUG: InitUserMeshData reached. build=%s %s\n", __DATE__, __TIME__);
    std::printf("### DEBUG: refinement=%s numlevel=%d derefine_count=%d\n",
                pin->GetString("mesh","refinement").c_str(),
                pin->GetInteger("mesh","numlevel"),
                pin->GetOrAddInteger("mesh","derefine_count",1));
    std::fflush(stdout);
  }

  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes SN ejecta.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // ambient medium everywhere
        phydro->u(IDN,k,j,i) = rho0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = p0/(gammagas-1.0);
        if (pscalars != nullptr) {
          pscalars->r(I_EJ,k,j,i) = 0.0;
        }
      }
    }
  }
}



//----------------------------------------------------------------------------------------
//! \fn void SNOuterX1()
//  \brief Sets boundary condition on downstream outflow boundary

void SNOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  bool applyDiode;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        // let gas outflow freely
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
//! \fn void SNInnerX1()
//  \brief Sets inner boundary where SN ejecta is injected

void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real rho, vel, pres;
  

  // physical time since explosion [s], never earlier than t_phys_floor
  Real t_phys = t_phys0 + time;
  t_phys = std::max(t_phys, t_phys_floor);

  // Radius of the inner boundary (use face coordinate for the boundary)
  Real r_in = pco->x1f(il);

    // Homologous velocity at the boundary: v = r / t
  vel = r_in / t_phys;
  if (vel > vmax) vel = vmax;


  // keep the original velocity law:
  // vel(t) = 1 / (t/Rej + 1/vmax)
  //vel = 1.0/(t_phys/Rej + 1.0/vmax);
  Real exponent = -(vel - v0_tail) / vc_tail;
  if (vel < v0_tail) exponent = 0.0;  // constant core at rho(v0,t)

  rho = rho0_tail * std::pow(t0_tail / t_phys, 3) * std::exp(exponent);

  pres = gasEntropy * std::pow(rho, gammagas);

  if ( Rej==0.0 || Eej==0.0 || Mej==0.0 || vmax==0.0 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in sn.cpp SNInnerX1" << std::endl
        << "ejecta parameters not set" << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

        prim(IDN,k,j,il-i) = rho;
        prim(IM1,k,j,il-i) = vel;
        prim(IM2,k,j,il-i) = 0.0;
        prim(IM3,k,j,il-i) = 0.0;
        prim(IEN,k,j,il-i) = pres;
        if (pmb->pscalars != nullptr) {
          pmb->pscalars->r(I_EJ,k,j,il-i) = 1.0;
        }

      }
    }
  }
}

int RefinementCondition(MeshBlock *pmb){
  static int printed = 0;
  if (!printed) {
    printed = 1;
    std::printf("### DEBUG: RefinementCondition entered. rank=%d gid=%d lid=%d pscalars=%p\n",
                Globals::my_rank, pmb->gid, pmb->lid, (void*)pmb->pscalars);
    std::fflush(stdout);
  }

  if (pmb->pscalars == nullptr) return 0;
  
  AthenaArray<Real> &u = pmb->phydro->u;
  AthenaArray<Real> &ej = pmb->pscalars->r;

  // current ejecta edge estimate
  Real time   = pmb->pmy_mesh->time;
  Real t_phys = t_phys0 + time;
  Real r_edge = vmax * t_phys;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        // only consider cells in a band around the ejecta front
        Real r_cell = pmb->pcoord->x1v(i);
        if (r_cell < 0.4 * r_edge || r_cell > 1.3 * r_edge) continue;

        if (i >= pmb->ie) continue;  // need a radial neighbor

        const Real ej_here = ej(I_EJ,k,j,i);
        const Real ej_out  = ej(I_EJ,k,j,i+1);
        const Real rho_here = u(IDN,k,j,i);
        const Real rho_out  = u(IDN,k,j,i+1);

      
         // --- criterion 1: tracer jump (ejecta edge) ---
        const Real dej = std::abs(ej_out - ej_here);
        if (dej > refine_dej_jump) return 1;

        // --- criterion 2: density cliff within ejecta (photosphere) ---
        if (ej_here > 0.5 && ej_out > 0.5) {
          const Real ratio = rho_here / rho_out;
          const Real r = (ratio > 1.0) ? ratio : 1.0 / ratio;
          if (r > refine_drho_frac) return 1;
        }
      }
    }
  }
  Real r_lo = pmb->pcoord->x1v(pmb->is);
  Real r_hi = pmb->pcoord->x1v(pmb->ie);
  if (r_hi < 0.3 * r_edge || r_lo > 1.5 * r_edge) return -1;
  
  return 0;  // no refinement needed
}



Real MyTimeStep(MeshBlock *pmb) {
  Real time   = pmb->pmy_mesh->time;
  Real dt     = pmb->pmy_mesh->dt;
  Real min_dt = largeDt;

  if ( time < smallDtDuration ) {
    min_dt = std::min( dt, smallDt );
  }

  return min_dt;
}
