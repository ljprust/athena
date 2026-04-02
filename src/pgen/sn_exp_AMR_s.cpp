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

//specific entropy parameters 
Real s0_fit, s1_fit;

// derived normalization A = rho0_tail * t0_tail^3
Real rho_norm_tail;

// minimum physical time (seconds) used in the fit
Real t_phys_floor;

// physical time at simulation t=0 (seconds)
Real t_phys0; 

constexpr int I_EJ = 0; //ejecta tracer

Real refine_drho_frac;
Real refine_dej_jump;

const Real a = 7.56e-15;
const Real mu = 4.0/3.0;
const Real mProton = 1.6726e-24;
const Real kB = 1.3807e-16;


static inline Real PressureFromTRho(Real temp, Real rho) {
 
  return (1.0/3.0)*a*std::pow(temp,4.0) + rho*kB*temp/(mu*mProton);
}
static inline Real EgasFromTRho(Real T, Real rho) {
  return a*std::pow(T, 4.0) + 1.5*rho*kB*T/(mu*mProton);
}

static inline Real EntropyFromTRho(Real T, Real rho) {
  const Real c = kB/(mu*mProton);
  return c*std::log(std::pow(T, 1.5)/rho) + (4.0/3.0)*a*std::pow(T, 3.0)/rho;
}

static Real TemperatureFromEntropyRho(Real s_target, Real rho) {
  auto f = [&](Real T) -> Real {
    return EntropyFromTRho(T, rho) - s_target;
  };

  // robust logarithmic bracket
  Real Tlo = 1.0e2;
  Real Thi = 1.0e8;

  while (f(Tlo) > 0.0) Tlo *= 0.5;
  while (f(Thi) < 0.0) Thi *= 2.0;

  // log-space bisection; monotonic because ds/dT > 0
  for (int n = 0; n < 80; ++n) {
    Real Tmid = std::sqrt(Tlo * Thi);
    if (f(Tmid) > 0.0) Thi = Tmid;
    else               Tlo = Tmid;
  }

  return std::sqrt(Tlo * Thi);
}

}

//takes in temperature and density and outputs pressure (from radIdeal EOS)


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
  t_phys_floor = pin->GetOrAddReal("problem","t_phys_floor",1.0);
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

  // new specific entropy fit parameters
  s0_fit = pin->GetOrAddReal("problem","s0_fit", 7.15e8);
  s1_fit = pin->GetOrAddReal("problem","s1_fit", 0.528);  

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
        //phydro->u(IEN,k,j,i) = p0/(gammagas-1.0);

        Real egas_gamma = p0/(gammagas-1.0);
        Real egas_eos   = peos->EgasFromRhoP(rho0, p0);

        phydro->u(IEN,k,j,i) = egas_eos;


        if (Globals::my_rank == 0 && i == is && j == js && k == ks) {
          std::printf("ambient: gamma_guess=%e eos=%e rel_diff=%e\n",
            egas_gamma, egas_eos,std::abs(egas_gamma-egas_eos)/std::max(egas_eos,1.0e-99));
          std::fflush(stdout);
}
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

  //pres = gasEntropy * std::pow(rho, gammagas);
  Real s_target = s0_fit + s1_fit * vel;
  const Real Tinj     = TemperatureFromEntropyRho(s_target, rho);
  const Real egas     = EgasFromTRho(Tinj, rho);
  pres  = pmb->peos->PresFromRhoEg(rho, egas);
  //pres = PressureFromEntropyRho(s_target, rho);

  //DEBUG: 
  const Real egas_chk = pmb->peos->EgasFromRhoP(rho, pres);
  const Real relerr = std::abs(egas_chk - egas) / std::max(egas, 1.0e-99);
  if (Globals::my_rank == 0) {
    std::printf("injection: t_phys=%e vel=%e rho=%e pres=%e egas=%e egas_chk=%e relerr=%e\n",
      t_phys, vel, rho, pres, egas, egas_chk, relerr);
    std::fflush(stdout);
  }

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

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real ej_here = ej(I_EJ,k,j,i);
        // check ejecta-side only
        if (ej_here <= 0.5) continue;
        const Real rho_here = u(IDN,k,j,i);    
        
        auto should_refine_with = [&](int kn, int jn, int in) -> bool {
          const Real ej_nb  = ej(I_EJ,kn,jn,in);
          const Real rho_nb = u(IDN,kn,jn,in);

          // 1) must look like an interface in tracer
          const Real dej = std::abs(ej_nb - ej_here);
          if (dej < refine_dej_jump) return false;

          // 2) must also have strong density contrast across that interface
          const Real diff = std::abs(rho_nb - rho_here);
          const Real ref  = std::max(rho_nb, rho_here);
          return diff > refine_drho_frac * ref;
        };

        if (i < pmb->ie && should_refine_with(k,j,i+1)) return 1;
        if (j < pmb->je && should_refine_with(k,j+1,i)) return 1;
        if (k < pmb->ke && should_refine_with(k+1,j,i)) return 1;
      }
    }
  }
  
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
