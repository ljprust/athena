// Model of companion star being hit by supernova

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // exp, pow, sqrt
#include <fstream>    // ifstream
#include <string>     // getline, string

// Athena++ headers
#include "../athena.hpp"                   // macros, enums, function signatures
#include "../athena_arrays.hpp"            // AthenaArray
#include "../globals.hpp"                  // Globals
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../hydro/hydro.hpp"              // Hydro
#include "../scalars/scalars.hpp"          // PassiveScalars
#include "../gravity/gravity.hpp"          // Gravity
#include "../mesh/mesh.hpp"                // mesh

// MPI header if needed
#ifdef MPI_PARALLEL
#include <mpi.h>  // MPI_COMM_WORLD
#endif

// OPENMP header if needed
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// Declarations
void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
		 const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
		 const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		 AthenaArray<Real> &cons_scalar);
void OutflowIX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
void OutflowOX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
void OutflowIX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
void OutflowOX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
void OutflowIX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
void OutflowOX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
		FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
		int ngh);
int GradientRefine(MeshBlock *pmb);
int RefinementCondition(MeshBlock *pmb);
Real MyTimeStep(MeshBlock *pmb);

// other functions
Real get_r(Real x, Real y, Real z, Real x0, Real y0);  // get distance from x0, y0 given (x,y,z)
Real get_rho_ejecta(Real v, Real t);           // get ejecta profile

// File declarations
namespace {
void DampVelocity(MeshBlock *pmb, const Real time, const Real dt,
    AthenaArray<Real> &cons);
void CreateEjecta(MeshBlock *pmb, const Real time, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar);
// void RotFrame(MeshBlock *pmb, const Real time, const Real dt,
// 	      const AthenaArray<Real> *flux,
// 	      const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
// 	      AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);
void ShiftFrame(MeshBlock *pmb, const Real time, AthenaArray<Real> &cons);
//void SourceMask(AthenaArray<Real> &src, int is, int ie, int js, int je,
//                int ks, int ke, const MGCoordinates &coord, Real time);
}

// File variables - problem definition
namespace {
Real rho_0;         // density scale (g/cm^3)
Real p_0;           // pressure scale (dyne/cm^2)

Real t_damp;     // time up to which velocity damping is in effect
Real tau_damp;   // velocity damping time scale
  
Real Prho53;        // P/rho**(5/3) from Ken Shen's data
Real a_sep;         // binary separation
Real v_orb;         // orbital velocity
Real v_orb_donor;   // orbital velocity of donor
Real v_m;           // maximum velocity (defines shock front)
Real r_m0;          // initial radius of shock front relative to explosion center
Real r_m;           // current radius of shock front relative to explosion center
Real t_ej;          // time after which ejecta is created
Real t_exp;         // explosion time
Real x_exp;         // explosion x-coordinate
Real t_since_exp;   // time after explosion
Real t_dur;         // ejecta source term duration

Real M_ej;          // total ejecta mass
Real E_ej;          // total explosion energy
int ejecta_case;    // which ejecta profile to use: exponential or gaussian
Real rho0_ej;       // normalization
Real v0_ej;         // velocity in Gaussian/exponential

Real M_d;           // donor mass
Real Omega_orb;     // orbital frequency
Real x_CMS;         // location of binary center of mass
Real t_roche_on;    // turn on Roche potential if time < t_roche_on [code]
  
Real phi_thres;  // threshold in grav potential in defining bulk motion
Real r0_thres;   // threshold in r0 in defining bulk motion
Real bulk_vx;    // vx of bulk
Real bulk_vy;    // vy of bulk
Real bulk_vz;    // vz of bulk
Real x_star;     // x of C.O.M. of star
Real y_star;     // y of C.O.M. of star
Real z_star;     // z of C.O.M. of star

Real t_freeze;   // code time to shift to COM frame
Real dt_freeze;  // time frame for shifting to COM frame

Real rho_fluff;  // density of background material in initial domain
Real p_fluff;    // pressure of background material in initial domain

Real r_mask;     // radius for gravity source mask
Real r_mask2;
Real t_mask;     // time to apply gravity source mask
  
Real t_initial;  // time before which apply a timestep limit
Real dt_initial; // apply this timestep initially

//Static mesh refinement -Sunny's attempt to combine the two 
//Real t_refine;   // time before which we apply an extra refine around donor
//Real x_refine;   // apply extra refinement for x in this range around origin
//Real y_refine;   // apply extra refinement for y/z in this range around origin
}

// File variables - physical unit scales
namespace {
Real x_0;    // length scale (cm)
Real t_0;    // time scale (s)
Real m_0;    // mass scale (g)
Real v_0;    // velocity scale (cm/s)
Real gg_0;   // gravitational constant scale (g.cm^3/s^2)
Real gg;     // gravitational constant in code units
}

// File variables - adaptive mesh refinement
namespace {
Real refine_grad_val;     // gradient threshold for refining
Real refine_grad_frac;    // gradient fraction for refining
Real derefine_grad_val;   // gradient threshold for derefining
Real derefine_grad_frac;  // gradient fraction for derefining
Real maxeps_threshold;    // gradient threshold for refining
Real r0_refine;           // donor concentration threshold for refining
Real r1_refine;           // ejecta concentration threshold for refining
Real rho_refine;          // density (code units) threshold for refining
}

// File variables - bookkeeping
namespace {
int num_rows_donor;
}

//----------------------------------------------------------------------------------------
// Function for preparing Mesh
// Inputs:
//   pin: input parameters
// Outputs: (none)

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Physical parameters
  const Real gg_cgs = 6.6743e-8;
  const Real rsun_cgs = 6.957e10;
  const Real kms_cgs = 1.0e5;
  const Real msun_cgs = 1.98840987e33;

  // Read scale parameters
  x_0 = pin->GetReal("problem", "x_0") * rsun_cgs;
  rho_0 = pin->GetReal("problem", "rho_0");
  p_0 = pin->GetReal("problem", "p_0");

  // Set physical scales
  v_0 = std::sqrt(p_0 / rho_0);
  t_0 = x_0 / v_0;
  m_0 = x_0 * SQR(x_0) * rho_0;
  gg_0 = p_0 / SQR(x_0 * rho_0);

  // Read remaining problem parameters
  t_damp = pin->GetReal("problem", "t_damp");
  tau_damp = pin->GetReal("problem", "tau_damp");
  
  Prho53 = pin->GetReal("problem", "Prho53")/p_0*std::pow(rho_0,5./3.);
  a_sep = pin->GetReal("problem", "a_sep") * rsun_cgs / x_0;
  v_orb = pin->GetReal("problem", "v_orb") * kms_cgs / v_0;

  t_ej = pin->GetReal("problem", "t_ej");
  t_dur = pin->GetReal("problem", "t_dur");
  
  M_ej = pin->GetReal("problem", "M_ej") * msun_cgs / m_0;
  E_ej = pin->GetReal("problem", "E_ej") / (m_0 * v_0 * v_0);
  ejecta_case = pin->GetInteger("problem", "ejecta_case");

  v_m = pin->GetReal("problem", "v_m") * kms_cgs / v_0;
  r_m0 = pin->GetReal("problem", "r_m0");

  phi_thres = pin->GetReal("problem", "phi_thres");
  r0_thres = pin->GetReal("problem", "r0_thres");

  // for shifting to COM frame
  t_freeze = pin->GetReal("problem", "t_freeze");
  dt_freeze = pin->GetReal("problem", "dt_freeze");

  // for gravity source mask
  r_mask = pin->GetReal("problem","r_mask");
  r_mask2 = r_mask * r_mask;
  t_mask = pin->GetReal("problem","t_mask");

  // for changing initial timestep
  t_initial = pin->GetReal("problem","t_initial");
  dt_initial = pin->GetReal("problem","dt_initial");

  // for ejecta profile
  rho0_ej = SQR(M_ej) * SQR(M_ej) * M_ej  / SQR(E_ej) / E_ej;

  // check ejecta case
  if (ejecta_case == 1) {
    // for Dwarkadas & Chevalier 1998 eqn 1
    rho0_ej *= SQR(6.0)*6.0;
    rho0_ej = std::sqrt(rho0_ej)/(8*PI);
    v0_ej = std::sqrt(E_ej/6.0/M_ej);
  } else {
    // Gaussian
    rho0_ej *= SQR(3.0/4.0/PI) * (3.0/4.0/PI);
    rho0_ej = std::sqrt(rho0_ej);
    v0_ej = std::sqrt(4.0/3.0*E_ej/M_ej);
  }


  // for Roche lobe
  M_d = pin->GetReal("problem", "M_d") * msun_cgs / m_0;

  // orbital velocity of donor
  v_orb_donor = v_orb * M_ej/(M_ej+M_d);

  t_roche_on = pin->GetReal("problem", "t_roche_on");

  // for fluff material
  rho_fluff = pin->GetReal("problem", "rho_fluff");
  p_fluff = pin->GetReal("problem", "p_fluff");
  
  // Read donor profile
  std::string data_file_donor = pin->GetString("problem", "data_file_donor");
  int num_comments_donor = pin->GetInteger("problem", "num_comments_donor");
  num_rows_donor = pin->GetInteger("problem", "num_rows_donor");
  int num_cols_donor = pin->GetInteger("problem", "num_cols_donor");
  int col_r_donor = pin->GetInteger("problem", "col_r_donor");
  int col_rho_donor = pin->GetInteger("problem", "col_rho_donor");
  int col_p_donor = pin->GetInteger("problem", "col_p_donor");

  // Read refinement criteria
  if (adaptive) {
    // refine_grad_val = pin->GetReal("problem", "refine_grad_val");
    // refine_grad_frac = pin->GetReal("problem", "refine_grad_frac");
    // derefine_grad_val = pin->GetReal("problem", "derefine_grad_val");
    // derefine_grad_frac = pin->GetReal("problem", "derefine_grad_frac");
    maxeps_threshold = pin->GetReal("problem", "maxeps_threshold");
    r0_refine = pin->GetReal("problem", "r0_refine");
    r1_refine = pin->GetReal("problem", "r1_refine");
    rho_refine = pin->GetReal("problem", "rho_refine");
    //t_refine = pin->GetReal("problem", "t_refine");
    //x_refine = pin->GetReal("problem", "x_refine");
    //y_refine = pin->GetReal("problem", "y_refine");
  }
  
  // Calculate explosion parameters
  t_exp = t_ej - r_m0 / v_m;
  x_exp = -a_sep;

  // Set gravitational parameters
  if (SELF_GRAVITY_ENABLED) {
    gg = gg_cgs / gg_0;
    SetGravitationalConstant(gg);
  }

  // Set parameters for rotating frame
  x_CMS = M_ej/(M_ej+M_d) * x_exp;
  Omega_orb = std::sqrt( gg*(M_ej+M_d)/a_sep/a_sep/a_sep );

  // Prepare arrays to hold profile
  AllocateRealUserMeshDataField(3);
  
  ruser_mesh_data[0].NewAthenaArray(num_rows_donor);
  ruser_mesh_data[1].NewAthenaArray(num_rows_donor);
  ruser_mesh_data[2].NewAthenaArray(num_rows_donor);

  // Alias arrays
  AthenaArray<Real> &r_star = ruser_mesh_data[0];
  AthenaArray<Real> &rho_star = ruser_mesh_data[1];
  AthenaArray<Real> &p_star = ruser_mesh_data[2];

  // Read profile data on rank 0
  if (Globals::my_rank == 0) {

    // Donor profile
    std::ifstream data_stream2(data_file_donor);
    std::string line;
    for (int n = 0; n < num_comments_donor; ++n) {
      std::getline(data_stream2, line);
    }
    Real vals2[num_cols_donor];
    for (int n = 0; n < num_rows_donor; ++n) {
      for (int m = 0; m < num_cols_donor; ++m) {
        data_stream2 >> vals2[m];
      }
      r_star(num_rows_donor - 1 - n) = vals2[col_r_donor] / x_0;
      rho_star(num_rows_donor - 1 - n) = vals2[col_rho_donor] / rho_0;
      p_star(num_rows_donor - 1 - n) = vals2[col_p_donor] / p_0;
    }
  }

  // Broadcast profile data to other ranks
  #ifdef MPI_PARALLEL
  {
    MPI_Bcast(r_star.data(), num_rows_donor, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho_star.data(), num_rows_donor, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(p_star.data(), num_rows_donor, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  }
  #endif

  // Enroll source term
  EnrollUserExplicitSourceFunction(SourceTerms);

  // Enroll boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, OutflowIX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, OutflowIX2);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x2, OutflowOX2);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, OutflowIX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, OutflowOX3);

  // Enroll time step
  EnrollUserTimeStepFunction(MyTimeStep);

  // Enroll gravity source term mask
  //EnrollUserMGGravitySourceMaskFunction(SourceMask);

  // Enroll refinement condition
  if (adaptive) {
    // EnrollUserRefinementCondition(GradientRefine);
    EnrollUserRefinementCondition(RefinementCondition);
  }

  
  return;
}

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters (not used)
// Outputs: (none)

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js - (ncells2 > 1 ? NGHOST : 0);
  int ju = je + (ncells2 > 1 ? NGHOST : 0);
  int kl = ks - (ncells3 > 1 ? NGHOST : 0);
  int ku = ke + (ncells3 > 1 ? NGHOST : 0);

  // Extract floors
  Real rho_floor = peos->GetDensityFloor();
  Real p_floor = peos->GetPressureFloor();

  // Alias arrays
  AthenaArray<Real> &r_star = pmy_mesh->ruser_mesh_data[0];
  AthenaArray<Real> &rho_star = pmy_mesh->ruser_mesh_data[1];
  AthenaArray<Real> &p_star = pmy_mesh->ruser_mesh_data[2];
  
  // Initialize primitive values
  for (int k = kl; k <= ku; ++k) {
    Real z = pcoord->x3v(k);
    for (int j = jl; j <= ju; ++j) {
      Real y = pcoord->x2v(j);

      Real x, r, rho_val, p_val;
      int index;
      // #ifdef OPENMP_PARALLEL
      // #pragma omp simd private(x,r,rho_val,p_val,index)
      // #endif
      for (int i = il; i <= iu; ++i) {

        // Determine radial location
        x = pcoord->x1v(i);
        r = std::sqrt(SQR(x) + SQR(y) + SQR(z));

        // Interpolate rho and p from tables
        rho_val = rho_fluff;
        p_val = p_fluff;
        if (r < r_star(num_rows_donor - 1)) {
          for (index = 1; index < num_rows_donor; ++index) {
            if (r < r_star(index)) {
              break;
            }
          }
          rho_val = rho_star(index - 1) + (rho_star(index) - rho_star(index - 1))
            / (r_star(index) - r_star(index - 1)) * (r - r_star(index - 1));
          p_val = p_star(index - 1) + (p_star(index) - p_star(index - 1))
            / (r_star(index) - r_star(index - 1)) * (r - r_star(index - 1));
        }
    
        // Set values
        if (rho_val > rho_fluff and p_val > p_fluff) {
          phydro->w(IDN,k,j,i) = rho_val;
          phydro->w(IPR,k,j,i) = p_val;
        } else {
          phydro->w(IDN,k,j,i) = rho_fluff;
          phydro->w(IPR,k,j,i) = p_fluff;
        }
        phydro->w(IVX,k,j,i) = 0.0;
        phydro->w(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = 0.0;
        if (NSCALARS >= 2) {
          pscalars->r(0,k,j,i) = (rho_val > rho_fluff and p_val > p_fluff) ? 1.0 : 0.0;
          pscalars->r(1,k,j,i) = 0.0;
        }
      }
    }
  }
  
  // Initialize conserved values
  AthenaArray<Real> b;
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);
  if (NSCALARS >= 2) {
    peos->PassiveScalarPrimitiveToConserved(pscalars->r, phydro->u, pscalars->s, pcoord,
        il, iu, jl, ju, kl, ku);
  }
  
  return;
}

//----------------------------------------------------------------------------------------
// Mask the density outside the donor star
/*namespace{
void SourceMask(AthenaArray<Real> &src, int is, int ie, int js, int je,
                int ks, int ke, const MGCoordinates &coord, Real time) {
  for (int k=ks; k<=ke; ++k) {
    Real z = coord.x3v(k);
    for (int j=js; j<=je; ++j) {
      Real y = coord.x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x = coord.x1v(i);
        Real r2 = x*x + y*y + z*z;
        if (r2 > r_mask2 and time < t_mask) {
          src(k, j, i) = 0.0;
        }
      }
    }
  }
  return;
}
}
*/
//----------------------------------------------------------------------------------------
Real get_rho_ejecta(Real v, Real t_since_exp) {
  // expoential form
  if (ejecta_case == 1) {
    return rho0_ej * std::exp(-v/v0_ej) / SQR(t_since_exp) / t_since_exp;
  } else {
    // Gaussian form
    return rho0_ej * std::exp(-SQR(v/v0_ej)) / SQR(t_since_exp) / t_since_exp;
  }
}

//----------------------------------------------------------------------------------------
// Collection of all source terms
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation
//   dt: timestep
//   prim: previous primitives (not used)
//   prim_scalar: previous primitive passive scalars (not used)
//   bcc: current magnetic field (not used)
//   cons: current conserved values
//   cons_scalar: current conserved passive scalars
// Outputs:
//   cons: updated with effect of sources
//   cons_scalar: updated with effect of sources

void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
		 const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
		 const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		 AthenaArray<Real> &cons_scalar) {

  // // Apply fictitious forces in rotating frame
  // if (time < t_roche_on) {
  //   RotFrame(pmb, time, dt, flux, prim, prim_scalar, cons, cons_scalar);
  //}

  // // Apply damping
  //if (time < t_damp) {
  //  DampVelocity(pmb, time, dt, cons);
  //}

  // // Create supernova ejecta
  if (time >= t_ej and time < t_ej + t_dur) {
    CreateEjecta(pmb, time, cons, cons_scalar);
  }

  // shift to COM frame (only for integrator stage 2 of vl2)
  /*if (stage == 2) {
    if (time >= t_freeze and time < t_freeze + dt_freeze) {
      ShiftFrame(pmb, time, cons);
    }
  }
  */
  
  return;
}

namespace{
  void ShiftFrame(MeshBlock *pmb, const Real time, AthenaArray<Real> &cons) {

    // Extract indices
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
	      Real rho = cons(IDN,k,j,i);
        Real ke_old = 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i))
            + SQR(cons(IM3,k,j,i))) / rho;
        Real ke_new = 0.5 * (SQR(cons(IM1,k,j,i)-rho*bulk_vx) + SQR(cons(IM2,k,j,i)-rho*bulk_vy)
            + SQR(cons(IM3,k,j,i)-rho*bulk_vz)) / rho;
        cons(IM1,k,j,i) += -rho*bulk_vx;
        cons(IM2,k,j,i) += -rho*bulk_vy;
        cons(IM3,k,j,i) += -rho*bulk_vz;
        cons(IEN,k,j,i) += ke_new - ke_old;
      }
    }
  }
  return;

  }
}

//----------------------------------------------------------------------------------------
// Velocity-damping source term
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation
//   dt: timestep
//   cons: current conserved values
// Outputs;
//   cons: updated with effect of sources
// Notes:
//   Only has an effect near the beginning of the simulation.

namespace {
void DampVelocity(MeshBlock *pmb, const Real time, const Real dt,
    AthenaArray<Real> &cons) {

  // Calculate damping factor
  Real dt_eff = std::exp(-time / (t_damp - time)) * dt;
  Real factor = std::exp(-dt_eff / tau_damp);

  // Extract indices
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real ke_old = 0.5 * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i))
            + SQR(cons(IM3,k,j,i))) / cons(IDN,k,j,i);
        Real ke_new = ke_old * SQR(factor);
        cons(IM1,k,j,i) *= factor;
        cons(IM2,k,j,i) *= factor;
        cons(IM3,k,j,i) *= factor;
        cons(IEN,k,j,i) += ke_new - ke_old;
      }
    }
  }
  return;
}
}

//----------------------------------------------------------------------------------------
// Co-rotating frame source term
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation
//   dt: timestep
//   cons: current conserved values
// Outputs;
//   cons: updated with effect of sources

namespace {
void RotFrame(MeshBlock *pmb, const Real time, const Real dt,
	      const AthenaArray<Real> *flux,
	      const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
	      AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {

  // Extract indices
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Go through cells
  for (int k = ks; k <= ke; ++k) {

    Real z = pmb->pcoord->x3v(k);
    Real zp1 = pmb->pcoord->x3v(k+1);
    Real zm1 = pmb->pcoord->x3v(k-1);

    Real dx3 = pmb->pcoord->dx3v(k);
    
    for (int j = js; j <= je; ++j) {

      Real y = pmb->pcoord->x2v(j);
      Real yp1 = pmb->pcoord->x2v(j+1);
      Real ym1 = pmb->pcoord->x2v(j-1);

      Real dx2 = pmb->pcoord->dx2v(j);
      
      for (int i = is; i <= ie; ++i) {
	
        Real x = pmb->pcoord->x1v(i);
        Real xp1 = pmb->pcoord->x1v(i+1);
        Real xm1 = pmb->pcoord->x1v(i-1);
        
        Real dx1 = pmb->pcoord->dx1v(i);
        
        // use donor material for rho
        Real rho = prim_scalar(0,k,j,i)*prim(IDN,k,j,i);
        // Real rho = prim(IDN,k,j,i);
        
        Real vx = prim(IVX,k,j,i);
        Real vy = prim(IVY,k,j,i);
        Real vz = prim(IVZ,k,j,i);

        // centrifugal force
        // Real a_cen_x = SQR(Omega_orb)*(x-x_CMS);
        // Real a_cen_y = SQR(Omega_orb)*y;
        // Real IEN_cen_x = prim(IDN,k,j,i)* vx * a_cen_x;
        // Real IEN_cen_y = prim(IDN,k,j,i)* vy * a_cen_y;
        
        Real OMG2 = SQR(Omega_orb);

        // for momentum term, average over cell volume
        Real a_cen_x = 0.5*(SQR(pmb->pcoord->x1f(i+1)) - SQR(pmb->pcoord->x1f(i)))/pmb->pcoord->dx1f(i) - x_CMS;
        Real a_cen_y = 0.5*(SQR(pmb->pcoord->x2f(j+1)) - SQR(pmb->pcoord->x2f(j)))/pmb->pcoord->dx2f(j);
        a_cen_x *= OMG2;
        a_cen_y *= OMG2;

        // for energy term, use eqn 57 of Mullen, Hanawa & Gammie
        Real phic = -0.5*OMG2*(SQR(x-x_CMS)+SQR(y));

        Real phil_x = 0.5*( -0.5*OMG2*(SQR(xm1-x_CMS)+SQR(y)) + phic );
        Real phir_x = 0.5*( -0.5*OMG2*(SQR(xp1-x_CMS)+SQR(y)) + phic );

        Real phil_y = 0.5*( -0.5*OMG2*(SQR(x-x_CMS)+SQR(ym1)) + phic );
        Real phir_y = 0.5*( -0.5*OMG2*(SQR(x-x_CMS)+SQR(yp1)) + phic );

        Real IEN_cen_x = -(flux[X1DIR](IDN,k,j,i)*(phic - phil_x) + flux[X1DIR](IDN,k,j,i+1)*(phir_x - phic))/dx1;
        Real IEN_cen_y = -(flux[X2DIR](IDN,k,j,i)*(phic - phil_y) + flux[X2DIR](IDN,k,j+1,i)*(phir_y - phic))/dx2;

        // coriolis force
        // Real psi = Omega_orb*dt;
        // Real a_cor_x = 2.0*Omega_orb*(-psi*vx + vy)/(1.0+SQR(psi));
        // Real a_cor_y = -2.0*Omega_orb*(vx + psi*vy)/(1.0+SQR(psi));

        // calculate accretor gravity
        Real GM = gg*M_ej;

        phic = -GM/get_r(x,y,z,x_exp,0.0);

        phil_x = 0.5*( -GM/get_r(xm1,y,z,x_exp,0.0) + phic );
        phir_x = 0.5*( -GM/get_r(xp1,y,z,x_exp,0.0) + phic );

        phil_y = 0.5*( -GM/get_r(x,ym1,z,x_exp,0.0) + phic );
        phir_y = 0.5*( -GM/get_r(x,yp1,z,x_exp,0.0) + phic );

        Real phil_z = 0.5*( -GM/get_r(x,y,zm1,x_exp,0.0) + phic );
        Real phir_z = 0.5*( -GM/get_r(x,y,zp1,x_exp,0.0) + phic );
        
        Real IM_g_x = -(phir_x-phil_x)/dx1;
        Real IM_g_y = -(phir_y-phil_y)/dx2;
        Real IM_g_z = -(phir_z-phil_z)/dx3;

        Real IEN_g_x = -(flux[X1DIR](IDN,k,j,i)*(phic - phil_x) + flux[X1DIR](IDN,k,j,i+1)*(phir_x - phic))/dx1;
        Real IEN_g_y = -(flux[X2DIR](IDN,k,j,i)*(phic - phil_y) + flux[X2DIR](IDN,k,j+1,i)*(phir_y - phic))/dx2;
        Real IEN_g_z = -(flux[X3DIR](IDN,k,j,i)*(phic - phil_z) + flux[X3DIR](IDN,k+1,j,i)*(phir_z - phic))/dx3;

        // update momentum
        Real acc_x = IM_g_x + a_cen_x;
        Real acc_y = IM_g_y + a_cen_y;
        Real acc_z = IM_g_z;
        
        cons(IM1,k,j,i) += dt*rho*acc_x;
        cons(IM2,k,j,i) += dt*rho*acc_y;
        cons(IM3,k,j,i) += dt*rho*acc_z;

        // energy
        // Real v_dot_a = vx*acc_x + vy*acc_y + vz*acc_z;
        // Real a_dot_a = SQR(acc_x) + SQR(acc_y) + SQR(acc_z);
        // cons(IEN,k,j,i) += rho*(v_dot_a*dt + 0.5*a_dot_a*SQR(dt));

        cons(IEN,k,j,i) += prim_scalar(0,k,j,i)*dt*(IEN_g_x + IEN_g_y + IEN_g_z);
        cons(IEN,k,j,i) += prim_scalar(0,k,j,i)*dt*(IEN_cen_x + IEN_cen_y);
        
      }
    }
  }
  return;
}
}

//----------------------------------------------------------------------------------------
// Time step
// Inputs:
//   pmb: pointer to MeshBlock
// Outputs:
//   returned value: minimum timestep

Real MyTimeStep(MeshBlock *pmb)
{
  Real time = pmb->pmy_mesh->time;
  Real dt = pmb->pmy_mesh->dt;
  Real next_time = time + dt;
  Real min_dt = 1e2;
  
  // start simulation with small timestep
  if (time == 0.0) {
    min_dt = 5e-4;
  }

  // apply small timestep
  if (time <= t_initial) {
    min_dt = std::min(min_dt,dt_initial);
  }

  // actual current time is next_time
  // assume proposed dt is similar to current dt
  // then if current time is before explosion initiates
  // and proposed next time is after explosion 
  // then severely limit timestep
  if (next_time <= t_ej and next_time + 1.2*dt >= t_ej) {
    min_dt = 5e-4;
  }

  if (next_time <= t_ej and next_time > t_ej - 1e-1) {
    min_dt = 5e-4;
  }
  
  return min_dt;
}

//----------------------------------------------------------------------------------------
// get distance relative to accretor position (located at x0)
Real get_r(Real x, Real y, Real z, Real x0, Real y0) {
  return std::sqrt( SQR(x-x0) + SQR(y-y0) + SQR(z) );
}

//----------------------------------------------------------------------------------------
// Ejecta source term
// Inputs:
//   pmb: pointer to MeshBlock
//   time: time of simulation
//   cons: current conserved values
//   cons_scalar: current conserved passive scalars
// Outputs;
//   cons: updated with effect of sources
//   cons_scalar: updated with effect of sources
// Notes:
//   Only has an effect for a brief duration.

namespace {
void CreateEjecta(MeshBlock *pmb, const Real time, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar) {

  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  // assume that when we inject ejecta, the donor is unbound immediately
  // and we move into that particular frame
  // the explosion center then moves with a negative y-velocity
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

  // Go through cells
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    Real z = pmb->pcoord->x3v(k);
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real y = pmb->pcoord->x2v(j);

      Real x, r, v, vx, vy, vz, rho_val, p_val, u, ke;
      // #ifdef OPENMP_PARALLEL
      // #pragma omp simd private(x, r, v, vx, vy, vz, rho_val, p_val, u, ke)
      // #endif
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        // Calculate position
        x = pmb->pcoord->x1v(i);
        r = get_r(x,y,z,x_exp,y_exp);

        if (r >= r_m0) {
          continue;
        }

        // Calculate velocity
        v = r / t_since_exp;
        vx = (x - x_exp) / r * v;
        vy = (y - y_exp) / r * v;
        vz = z / r * v;
        
        rho_val = get_rho_ejecta(v,t_since_exp);

        if (rho_val <= rho_floor) {
          continue;
        }

        // Calculate energy
        // assume entropy ~ P/rho**(5/3) is constant
        // Get P from that
        if (GENERAL_EOS) {
          p_val = Prho53*std::pow(rho_val,5./3.);
          u = pmb->peos->EgasFromRhoP(rho_val, p_val);
        } else {
          u = 3./2.*Prho53*std::pow(rho_val,5./3.);
        }

        // Calculate effect of orbital velocity
        vy -= v_orb;
        ke = 0.5 * rho_val * (SQR(vx) + SQR(vy) + SQR(vz));

        // Set conserved variables
        cons(IDN,k,j,i) = rho_val;
        cons(IM1,k,j,i) = rho_val * vx;
        cons(IM2,k,j,i) = rho_val * vy;
        cons(IM3,k,j,i) = rho_val * vz;
        cons(IEN,k,j,i) = ke + u;
        if (NSCALARS >= 2) {
          cons_scalar(0,k,j,i) = 0.0;
          cons_scalar(1,k,j,i) = rho_val;
        }
      }
    }
  }
  return;
}
}



// refinement condition: check the maximum pressure gradient
// and maximum scalar gradient.
// refine only if donor scalar is greater than r0_refine. 
int RefinementCondition(MeshBlock *pmb) {

  // Prepare index bounds
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // coordinates
  Coordinates *pcoord = pmb->pcoord;

  // primitives
  AthenaArray<Real> &w = pmb->phydro->w;

  // Maximum gradient
  Real maxeps = 0.0;

  // current time
  Real time = pmb->pmy_mesh->time;

  // level of mesh block
  //int pmb_level = pmb->loc.level - pmb->pmy_mesh->root_level;

  // refine 1 level after taking 1 step
  // if (pmb_level == 0 and time > 1e-5) {
  //   return 1;
  // }

  // further refine using gradients
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {

        Real eps;
        if (NSCALARS>=2) {
          AthenaArray<Real> &r = pmb->pscalars->r;
          // refine if donor concentration exceeds threshold
          // if ((r(0,k,j,i) >= r0_refine) ||
          // (r(1,k,j,i) >= r1_refine and pmb->loc.level == pmb->pmy_mesh->root_level)) {
          // if (r(0,k,j,i) >= r0_refine) {
          if (r(0,k,j,i) >= r0_refine and w(IDN,k,j,i) >= rho_refine) {
            Real epsp = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                      +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                      +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
            Real epsr = std::sqrt(SQR(0.5*(w(IDN,k,j,i+1) - w(IDN,k,j,i-1)))
                      +SQR(0.5*(w(IDN,k,j+1,i) - w(IDN,k,j-1,i)))
                      +SQR(0.5*(w(IDN,k+1,j,i) - w(IDN,k-1,j,i))))/w(IDN,k,j,i);
            eps = std::max(epsp,epsr);
            maxeps = std::max(maxeps, eps);
          }


        } else {
          Real epsp = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                    +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                    +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
          Real epsr = std::sqrt(SQR(0.5*(w(IDN,k,j,i+1) - w(IDN,k,j,i-1)))
                    +SQR(0.5*(w(IDN,k,j+1,i) - w(IDN,k,j-1,i)))
                    +SQR(0.5*(w(IDN,k+1,j,i) - w(IDN,k-1,j,i))))/w(IDN,k,j,i);
          eps = std::max(epsp,epsr);
          maxeps = std::max(maxeps, eps);
        }
      }
    }
  }
  
  // refine some volume around the donor if time < t_refine
  /*
  if (pmb->pmy_mesh->time < t_refine) {
    if ((pcoord->x1f(is) >= -x_refine && pcoord->x1f(ie+1) <= x_refine) ||
    (pcoord->x1f(is) <= -x_refine && pcoord->x1f(ie+1) >= -x_refine) ||
    (pcoord->x1f(is) <= x_refine && pcoord->x1f(ie+1) >= x_refine) ) {
      if ((pcoord->x2f(js) >= -y_refine && pcoord->x2f(je+1) <= y_refine) ||
      (pcoord->x2f(js) <= -y_refine && pcoord->x2f(je+1) >= -y_refine) ||
      (pcoord->x2f(js) <= y_refine && pcoord->x2f(je+1) >= y_refine) ) {
        if ((pcoord->x3f(ks) >= -y_refine && pcoord->x3f(ke+1) <= y_refine) ||
        (pcoord->x3f(ks) <= -y_refine && pcoord->x3f(ke+1) >= -y_refine) ||
        (pcoord->x3f(ks) <= y_refine && pcoord->x3f(ke+1) >= y_refine) ) {
          // only refine if at root level
          if (pmb_level == 0) {
            return 1;
          }
        }
      }
    }
  } */

  

  if (maxeps > maxeps_threshold) return 1;
  if (maxeps < 0.25*maxeps_threshold) return -1;
  // if (pmb_level > 1) {
  //   if (maxeps < 0.25*maxeps_threshold) return -1;
  // }

  return 0;
}

//----------------------------------------------------------------------------------------
// Function called once every time step for user-defined work.
// Inputs: (none)
// Outputs: (none)
/*void Mesh::UserWorkInLoop() {
  
  return;
}
*/
//----------------------------------------------------------------------------------------
// Outflow boundary condition - left
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowIX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias arrays
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is - ngh; i <= is - 1; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,k,j,is);
        prim(IPR,k,j,i) = prim(IPR,k,j,is);
        prim(IVX,k,j,i) = std::min(prim(IVX,k,j,is), 0.0);
        prim(IVY,k,j,i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,i) = prim(IVZ,k,j,is);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
	
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

        Real rho_val = get_rho_ejecta(v,t_since_exp);

        if (r < r_m and rho_val > rho_floor and vx > 0.0) {

          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {

            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition - right
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowOX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias array of positions
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = ie + 1; i <= ie + ngh; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,k,j,ie);
        prim(IPR,k,j,i) = prim(IPR,k,j,ie);
        prim(IVX,k,j,i) = std::max(prim(IVX,k,j,ie), 0.0);
        prim(IVY,k,j,i) = prim(IVY,k,j,ie);
        prim(IVZ,k,j,i) = prim(IVZ,k,j,ie);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
	
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

        Real rho_val = get_rho_ejecta(v,t_since_exp);
        
        if (r < r_m and rho_val > rho_floor and vx < 0.0) {
          
          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {

            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition - front
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowIX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias array of positions
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js - ngh; j <= js - 1; ++j) {
      for (int i = is; i <= ie; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,k,js,i);
        prim(IPR,k,j,i) = prim(IPR,k,js,i);
        prim(IVX,k,j,i) = prim(IVX,k,js,i);
        prim(IVY,k,j,i) = std::min(prim(IVY,k,js,i), 0.0);
        prim(IVZ,k,j,i) = prim(IVZ,k,js,i);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
	
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

        Real rho_val = get_rho_ejecta(v,t_since_exp);
        
        if (r < r_m and rho_val > rho_floor and vy > 0.0) {
          
          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {
            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition - back
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowOX2(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias array of positions
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;

  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = je + 1; j <= je + ngh; ++j) {
      for (int i = is; i <= ie; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,k,je,i);
        prim(IPR,k,j,i) = prim(IPR,k,je,i);
        prim(IVX,k,j,i) = prim(IVX,k,je,i);
        prim(IVY,k,j,i) = std::max(prim(IVY,k,je,i), 0.0);
        prim(IVZ,k,j,i) = prim(IVZ,k,je,i);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
        
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

        Real rho_val = get_rho_ejecta(v,t_since_exp);
        
        if (r < r_m and rho_val > rho_floor and vy < 0.0) {
          
          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {

            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition - bottom
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowIX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias array of positions
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;

  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ks - ngh; k <= ks - 1; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,ks,j,i);
        prim(IPR,k,j,i) = prim(IPR,ks,j,i);
        prim(IVX,k,j,i) = prim(IVX,ks,j,i);
        prim(IVY,k,j,i) = prim(IVY,ks,j,i);
        prim(IVZ,k,j,i) = std::min(prim(IVZ,ks,j,i), 0.0);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
	
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

        Real rho_val = get_rho_ejecta(v,t_since_exp);

        if (r < r_m and rho_val > rho_floor and vz > 0.0) {
          
          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {

            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Outflow boundary condition - top
// Inputs:
//   pmb: pointer to MeshBlock
//   pcoord: pointer to Coordinates
//   prim: primitives
//   b: magnetic field (not used)
//   time: simulation time
//   dt: effective time step (not used)
//   is, ie, js, je, ks, ke: indices demarcating active region
//   ngh: effective number of ghost zones
// Outputs:
//   prim: primitives set in ghost zones
//   b: magnetic field (not set)
// Notes:
//   Primitive variables are extrapolated as constant but with no inflow most of the time;
//       otherwise they are set to SN ejecta in certain time range.

void OutflowOX3(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
    FaceField &bb, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,
    int ngh) {

  // Alias array of positions
  AthenaArray<Real> &x = pcoord->x1v;
  AthenaArray<Real> &y = pcoord->x2v;
  AthenaArray<Real> &z = pcoord->x3v;
  // AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  AthenaArray<Real> &prim_scalar = (prim.data() == pmb->phydro->coarse_prim_.data()) ? pmb->pscalars->coarse_r_ : pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real t_since_ej = time - t_ej;
  Real r_m = v_m * t_since_exp;
  Real y_exp = - v_orb * t_since_ej ; // use t_ej instead of t_exp
  
  Real rho_floor = pmb->peos->GetDensityFloor();

  // Go through cells
  for (int k = ke + 1; k <= ke + ngh; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {

        // Outflow as base case
        prim(IDN,k,j,i) = prim(IDN,ke,j,i);
        prim(IPR,k,j,i) = prim(IPR,ke,j,i);
        prim(IVX,k,j,i) = prim(IVX,ke,j,i);
        prim(IVY,k,j,i) = prim(IVY,ke,j,i);
        prim(IVZ,k,j,i) = std::max(prim(IVZ,ke,j,i), 0.0);
        if (NSCALARS >= 2) {
          prim_scalar(0,k,j,i) = 0.0;
          prim_scalar(1,k,j,i) = 0.0;
        }
	
        // Calculate position
        Real r = get_r(x(i),y(j),z(k),x_exp,y_exp);

        // Calculate velocity
        Real v = r / t_since_exp;
        Real vx = (x(i) - x_exp) / r * v;
        Real vy = (y(j) - y_exp) / r * v;
        Real vz = z(k) / r * v;

        // Calculate effect of orbital velocity
        vy -= v_orb;

	      Real rho_val = get_rho_ejecta(v,t_since_exp);

        if (r < r_m and rho_val > rho_floor and vz < 0.0) {
          
          // Calculate pressure
          // assume entropy ~ P/rho**(5/3) is constant
          // Get P from that
          Real p_val = Prho53*std::pow(rho_val,5./3.);
          Real p_floor = pmb->peos->GetPressureFloor();
          if (p_val > p_floor) {

            // Set primitive variables
            prim(IDN,k,j,i) = rho_val;
            prim(IPR,k,j,i) = p_val;
            prim(IVX,k,j,i) = vx;
            prim(IVY,k,j,i) = vy;
            prim(IVZ,k,j,i) = vz;
            if (NSCALARS >= 2) {
              prim_scalar(1,k,j,i) = 1.0;
            }
          }
        }
      }
    }
  }
  return;
}