//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sn_exp_AMR_s_table2.cpp 
//! \brief Initializes SN ejecta into polar wedge
//4/22/26
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
#include <vector>

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
Real v_c, rho_0;

//specific entropy parameters 
Real s0_fit, s1_fit;

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

//table Helmholz EOS stuff
bool use_temperature_table;
std::string temperature_table_file;

std::vector<Real> inj_t_table;
std::vector<Real> inj_T_table;

static void LoadTemperatureTable(const std::string& filename) {
  //Load temperature table from file. The file should have three columns: time [s], velocity [cm/s], temperature [K]
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    std::stringstream msg;
    msg << "### FATAL ERROR\n"
        << "Could not open temperature table file: " << filename << "\n";
    ATHENA_ERROR(msg);
  }

  inj_t_table.clear();
  inj_T_table.clear();

  Real tval, vval, Tval;
  while (fin >> tval >> vval >> Tval) {
    inj_t_table.push_back(tval);
    inj_T_table.push_back(Tval);
  }

  if (inj_t_table.size() < 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in sn.cpp\n"
        << "Temperature table must contain at least 2 rows.\n";
    ATHENA_ERROR(msg);
  }

  for (std::size_t i = 1; i < inj_t_table.size(); ++i) {
    if (inj_t_table[i] <= inj_t_table[i-1]) {
      std::stringstream msg;
      msg << "### FATAL ERROR in sn.cpp\n"
          << "Temperature table times must be strictly increasing.\n";
      ATHENA_ERROR(msg);
    }
  }

  if (Globals::my_rank == 0) {
    std::printf("Loaded temperature table %s with %zu rows\n",
                filename.c_str(), inj_t_table.size());
    std::printf("  t range: [%e, %e]\n", inj_t_table.front(), inj_t_table.back());
    std::fflush(stdout);
  }
}

static Real InterpolateTemperature(Real t) {
  // interpolator for the temperature table 
  if (inj_t_table.empty()) {
    std::stringstream msg;
    msg << "### FATAL ERROR \n"
        << "Temperature table requested but not loaded.\n";
    ATHENA_ERROR(msg);
  }

  // Clamp outside the table range
  if (t <= inj_t_table.front()) return inj_T_table.front();
  if (t >= inj_t_table.back())  return inj_T_table.back();

  auto it = std::upper_bound(inj_t_table.begin(), inj_t_table.end(), t);
  std::size_t i1 = static_cast<std::size_t>(it - inj_t_table.begin());
  std::size_t i0 = i1 - 1;

  Real t0 = inj_t_table[i0];
  Real t1 = inj_t_table[i1];
  Real T0 = inj_T_table[i0];
  Real T1 = inj_T_table[i1];

  Real w = (t - t0) / (t1 - t0);
  return (1.0 - w) * T0 + w * T1;
}



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

  //  tail-fit parameters (v_c, rho_0))
  v_c        = pin->GetOrAddReal("problem","v_c",1.2e8);
  rho_0      = pin->GetOrAddReal("problem","rho_0", 3.49e7);

  // earliest physical time represented in the fit (seconds)
  t_phys_floor = pin->GetOrAddReal("problem","t_phys_floor",1.0);
  t_phys0 = pin->GetOrAddReal("problem","t_phys0", Rej / vmax);

  // properties of initial timestep
  smallDt = pin->GetReal("problem", "smallDt");
  largeDt = pin->GetReal("problem", "largeDt");
  smallDtDuration = pin->GetReal("problem", "smallDtDuration");

  //refinement parameters
  refine_drho_frac = pin->GetOrAddReal("problem","refine_drho_frac",0.1);
  refine_dej_jump  = pin->GetOrAddReal("problem","refine_dej_jump",0.4);

  // new specific entropy fit parameters
  s0_fit = pin->GetOrAddReal("problem","s0_fit", 3.63e7);
  s1_fit = pin->GetOrAddReal("problem","s1_fit", 0.531);  

  // tell athena to use our user-defined boundaries
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);

  // set a very small timestep initially
  EnrollUserTimeStepFunction(MyTimeStep);

  //enroll refinement condition
  EnrollUserRefinementCondition(RefinementCondition);

  use_temperature_table =
    pin->GetOrAddBoolean("problem", "use_temperature_table", false);
  temperature_table_file =
      pin->GetOrAddString("problem", "temperature_table_file", "");

  if (use_temperature_table) {
    if (temperature_table_file.empty()) {
      std::stringstream msg;
      msg << "### FATAL ERROR\n"
          << "use_temperature_table=true but temperature_table_file is empty.\n";
      ATHENA_ERROR(msg);
    }
    LoadTemperatureTable(temperature_table_file);
  }


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
  Real exponent = -(vel) / v_c;
  //if (vel < v0_tail) exponent = 0.0;  // constant core at rho(v0,t)

  rho = rho_0 * std::pow(1.0 / t_phys, 3) * std::exp(exponent);
  //rho = std::max(rho, 1e-11); //atempting to clamp rho

  Real Tinj;
  if (use_temperature_table) {
    Tinj = InterpolateTemperature(time); //test to see if multiplication fixes numerical issue
  } else {
    Real s_target = s0_fit + s1_fit * vel;
    Tinj = TemperatureFromEntropyRho(s_target, rho);
  }
  //const Real egas   = EgasFromTRho(Tinj, rho); //broek because density was too small to invert
  //pres  = pmb->peos->PresFromRhoEg(rho, egas);
  //pres = PressureFromEntropyRho(s_target, rho);
  pres = PressureFromTRho(Tinj, rho);


    //const Real s_target = s0_fit + s1_fit * vel;
    //const Real s_inj_athena = EntropyFromTRho(Tinj, rho);
  const Real egas_direct = EgasFromTRho(Tinj, rho);
  const Real egas_eos = pmb->peos->EgasFromRhoP(rho, pres);

    if (Globals::my_rank == 0) {
      std::printf(
        "inj check: t_phys=%e vel=%e rho=%e Tinj=%e pres=%e "
        "s_target=%e s_inj_athena=%e s_ratio=%e "
        "egas_direct=%e egas_eos=%e e_ratio=%e\n",
        t_phys, vel, rho, Tinj, pres,
        egas_direct, egas_eos, egas_eos/std::max(egas_direct,1.0e-99)
      );
      std::fflush(stdout);
    }

  //DEBUG: 
  //const Real egas_chk = pmb->peos->EgasFromRhoP(rho, pres);
  //const Real relerr = std::abs(egas_chk - egas) / std::max(egas, 1.0e-99);
    const Real eint = EgasFromTRho(Tinj, rho);   // internal energy density
    const Real rhov2 = rho * vel * vel;          // here vel2=vel3=0 at injection
    const Real ratio = eint / std::max(rhov2, 1.0e-99);
    const Real pfloor   = pmb->peos->GetPressureFloor();


    const Real pres_direct = PressureFromTRho(Tinj, rho);
    const Real pres_rt     = pmb->peos->PresFromRhoEg(rho, egas_direct);
    const Real egas_rt     = pmb->peos->EgasFromRhoP(rho, pres_direct);


    

    if (Globals::my_rank == 0) {
      std::printf(
        "inj diag: t_phys=%e vel=%e rho=%e Tinj=%e pres=%e "
        "eint=%e rhov2=%e eint_over_rhov2=%e\n",
        t_phys, vel, rho, Tinj, pres, eint, rhov2, ratio
      );
      std::printf(
          "EOS RT: T=%e rho=%e p_dir=%e p_rt=%e p_ratio=%e p_floor= %e e_dir=%e e_rt=%e e_ratio=%e\n",
          Tinj, rho,
          pres_direct, pres_rt, pres_rt/std::max(pres_direct,1.0e-99), pfloor,
          egas_direct, egas_rt, egas_rt/std::max(egas_direct,1.0e-99));
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

          if (Globals::my_rank == 0 && j == jl && k == kl && i == 1) {
                std::printf("ghost : rho=%e v=%e p=%e\n",
                            prim(IDN,k,j,il-1), prim(IM1,k,j,il-1), prim(IEN,k,j,il-1));
                std::printf("active: rho=%e v=%e p=%e\n",
                            prim(IDN,k,j,il  ), prim(IM1,k,j,il  ), prim(IEN,k,j,il  ));
                std::fflush(stdout);
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

  // Pin the inner-X1 boundary block at root level. Prolongation across the
  // steep injection contact produces fine cells where C2P puts rho at dfloor
  // with NaN v,p; SNInnerX1 only writes ghost zones so the corrupted active
  // cell never heals. Refinement is only needed at the ejecta edge.
  if (pmb->loc.lx1 == 0) return -1;

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
