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

// Declarations
void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
		 const AthenaArray<Real> *flux,
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



// for history

Real GetDonorCOM(MeshBlock *pmb, int iout);   // get center of mass of donor
Real GetBulkV(MeshBlock *pmb, int iout);      // get velocity of star used to define Be

Real GetBoundMass(MeshBlock *pmb, int iout);  // get total bound mass (Be<0) (unlabeled, r0, r1)

Real GetBoundMom(MeshBlock *pmb, int iout);   // get momentum of bound material (Be<0)

Real GetDonorVol(MeshBlock *pmb, int iout);   // get volume of Be<0 material for bound remnant/donor
Real GetDonorR(MeshBlock *pmb, int iout);     // get min radius of Be<0 material for bound remnant/donor

Real GetDonorPDc(MeshBlock *pmb, int iout);   // central pressure/density in donor
Real GetMaxPD(MeshBlock *pmb, int iout);      // max pressure/density

Real GetM_v(MeshBlock *pmb, int iout);        // mass in certain speed bin/ theta bin

Real GetBoundMassPureSN(MeshBlock *pmb, int iout);  // get SN bound mass (Be<0) that is pure (r1 > 0.95)

Real GetMdotBound(MeshBlock *pmb, int iout);  // get donor bound mass flow (Be<0) through the box

Real GetMFloor(MeshBlock *pmb, int iout);     // get total mass that hits the density floor


// other functions
Real get_r(Real x, Real y, Real z, Real x0);  // get distance from x0 given (x,y,z)
Real get_rho_gauss(Real v, Real t);           // get ejecta profile

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
// void ShiftFrame(MeshBlock *pmb, const Real time, AthenaArray<Real> &cons);
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
  r_m0 = pin->GetReal("problem", "r_m0") * rsun_cgs / x_0;

  phi_thres = pin->GetReal("problem", "phi_thres");
  r0_thres = pin->GetReal("problem", "r0_thres");

  // for shifting to COM frame
  // t_freeze = pin->GetReal("problem", "t_freeze");

  // for ejecta profile
  // (3/4/PI)^{3/2} * M^5/2 / E^3/2 
  rho0_ej = SQR(M_ej) * SQR(M_ej) * M_ej  / SQR(E_ej) / E_ej;
  rho0_ej *= SQR(3.0/4.0/PI) * (3.0/4.0/PI);
  rho0_ej = std::sqrt(rho0_ej);

  v0_ej = std::sqrt(4.0/3.0*E_ej/M_ej);

  // for Dwarkadas & Chevalier 1998 eqn 1
  // rho0_ej = SQR(M_ej) * SQR(M_ej) * M_ej  / SQR(E_ej) / E_ej;
  // rho0_ej *= SQR(6.0)*6.0;
  // rho0_ej = std::sqrt(rho0_ej)/(8*PI);

  // v0_ej = std::sqrt(E_ej/6.0/M_ej);

  // for Roche lobe
  M_d = pin->GetReal("problem", "M_d") * msun_cgs / m_0;

  t_roche_on = pin->GetReal("problem", "t_roche_on");
  
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
    refine_grad_val = pin->GetReal("problem", "refine_grad_val");
    refine_grad_frac = pin->GetReal("problem", "refine_grad_frac");
    derefine_grad_val = pin->GetReal("problem", "derefine_grad_val");
    derefine_grad_frac = pin->GetReal("problem", "derefine_grad_frac");
    maxeps_threshold = pin->GetReal("problem", "maxeps_threshold");
    r0_refine = pin->GetReal("problem", "r0_refine");
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

  // Enroll refinement condition
  if (adaptive) {
    // EnrollUserRefinementCondition(GradientRefine);
    EnrollUserRefinementCondition(RefinementCondition);
  }

  // Enroll history output
  //AllocateUserHistoryOutput(201);

  // get center of mass of donor
 // EnrollUserHistoryOutput(0, GetDonorCOM, "DonorX", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(1, GetDonorCOM, "DonorY", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(2, GetDonorCOM, "DonorZ", UserHistoryOperation::max);

  // get velocity of star used to define Be
  //EnrollUserHistoryOutput(3, GetBulkV, "BulkVx", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(4, GetBulkV, "BulkVy", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(5, GetBulkV, "BulkVz", UserHistoryOperation::max);

  // get total bound mass (Be<0) (unlabeled, r0, r1)
  //EnrollUserHistoryOutput(6, GetBoundMass, "BoundMass", UserHistoryOperation::sum);
  //EnrollUserHistoryOutput(7, GetBoundMass, "BoundMass0", UserHistoryOperation::sum);
  //EnrollUserHistoryOutput(8, GetBoundMass, "BoundMass1", UserHistoryOperation::sum);

  // get momentum of bound material (Be<0)
  //EnrollUserHistoryOutput(9, GetBoundMom, "BoundMomx", UserHistoryOperation::sum);
  //EnrollUserHistoryOutput(10, GetBoundMom, "BoundMomy", UserHistoryOperation::sum);
  //EnrollUserHistoryOutput(11, GetBoundMom, "BoundMomz", UserHistoryOperation::sum);

  // get central pressure/density at the donor C.O.M.
  //EnrollUserHistoryOutput(12, GetDonorPDc, "DonorRhoc", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(13, GetDonorPDc, "DonorPc", UserHistoryOperation::max);

  // get max pressure/density
  //EnrollUserHistoryOutput(14, GetMaxPD, "DonorMaxRho", UserHistoryOperation::max);
  //EnrollUserHistoryOutput(15, GetMaxPD, "DonorMaxP", UserHistoryOperation::max);

  // get outflowing mass in certain velocity range

  //Real v_bin = 9.30;
  //for (int i = 0; i <= 46; i++) {
    //std::string num_text = std::to_string(v_bin);
    //std::string rounded = num_text.substr(0, num_text.find(".")+3);
    //rounded.replace(rounded.find("."),1,"p");
    // std::cout << rounded << std::endl;
    //std::string myhead;
    //myhead = "v0_";
    //myhead.append(rounded);

    //EnrollUserHistoryOutput(16+i, GetM_v, myhead.c_str(), UserHistoryOperation::sum);

    //v_bin -= 0.05;
  //}

  //v_bin = 9.30;
  //for (int i = 0; i <= 46; i++) {
    //std::string num_text = std::to_string(v_bin);
    //std::string rounded = num_text.substr(0, num_text.find(".")+3);
    //rounded.replace(rounded.find("."),1,"p");
    // std::cout << rounded << std::endl;
    //std:: string myhead = "v1_";
    //myhead.append(rounded);

    //EnrollUserHistoryOutput(63+i, GetM_v, myhead.c_str(), UserHistoryOperation::sum);

    //v_bin -= 0.05;
 // }

  //EnrollUserHistoryOutput(110, GetMdotBound, "MdotBound", UserHistoryOperation::sum);

  // get outflowing mass in certain theta range
  // for donor material
  //Real theta_bin = 0.0;
  //for (int i = 0; i < 45; i++) {
    //std::string num_text = std::to_string(theta_bin);
    //std::string rounded = num_text.substr(0, num_text.find(".")+3);
    //rounded.replace(rounded.find("."),1,"p");
        
    //std::string myhead = "theta0_";
    //myhead.append(rounded);
    
    // std::cout << myhead << std::endl;
    //EnrollUserHistoryOutput(111+i, GetM_v, myhead.c_str(), UserHistoryOperation::sum);

    //theta_bin += 2.0;
  //}

  // now for SN ejecta
  //theta_bin = 0.0;
  //for (int i = 0; i < 45; i++) {
    //std::string num_text = std::to_string(theta_bin);
    //std::string rounded = num_text.substr(0, num_text.find(".")+3);
    //rounded.replace(rounded.find("."),1,"p");
        
    //std::string myhead = "theta1_";
    //myhead.append(rounded);
    
    // std::cout << myhead << std::endl;
    //EnrollUserHistoryOutput(156+i, GetM_v, myhead.c_str(), UserHistoryOperation::sum);

    //theta_bin += 2.0;
  //}
  
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
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {

	// Determine radial location
	Real x = pcoord->x1v(i);
	Real y = pcoord->x2v(j);
	Real z = pcoord->x3v(k);
	
	Real r = std::sqrt(SQR(x) + SQR(y) + SQR(z));
	    

	// Interpolate rho and p from tables
	Real rho_val = rho_floor;
	Real p_val = p_floor;
	if (r < r_star(num_rows_donor - 1)) {
	  int index;
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
	
	rho_val = std::max(rho_val, rho_floor);
	p_val = std::max(p_val, p_floor);
	
	// Set values
	phydro->w(IDN,k,j,i) = rho_val;
	phydro->w(IPR,k,j,i) = p_val;
	phydro->w(IVX,k,j,i) = 0.0;
	phydro->w(IVY,k,j,i) = 0.0;
	phydro->w(IVZ,k,j,i) = 0.0;
	if (NSCALARS >= 2) {
	  pscalars->r(0,k,j,i) = (rho_val > rho_floor and p_val > p_floor) ? 1.0 : 0.0;
	  pscalars->r(1,k,j,i) = 0.0;
	}
      }
    }
  }
  
    
  
  // Initialize conserved values
  AthenaArray<Real> b;
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);
  if (NSCALARS > 0) {
    peos->PassiveScalarPrimitiveToConserved(pscalars->r, phydro->u, pscalars->s, pcoord,
        il, iu, jl, ju, kl, ku);
  }

  
  return;
}

//----------------------------------------------------------------------------------------
Real get_rho_gauss(Real v, Real t_since_exp) {

  return rho0_ej * std::exp(-SQR(v/v0_ej)) / SQR(t_since_exp) / t_since_exp;
  
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
		 const AthenaArray<Real> *flux,
		 const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
		 const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
		 AthenaArray<Real> &cons_scalar) {

  // // Apply fictitious forces in rotating frame
  // if (time < t_roche_on) {
  //   RotFrame(pmb, time, dt, flux, prim, prim_scalar, cons, cons_scalar);
  // }

  // // Apply damping
  if (time < t_damp) {
    DampVelocity(pmb, time, dt, cons);
  }

  // // Create supernova ejecta
  if (time >= t_ej and time < t_ej + t_dur) {
    CreateEjecta(pmb, time, cons, cons_scalar);
  }

  // // shift to COM frame
  // if (time >= t_freeze and time < t_freeze + t_dur) {
  //   ShiftFrame(pmb, time, cons);
  // }
  
  return;
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
Real get_r(Real x, Real y, Real z, Real x0)
{
  return std::sqrt( SQR(x-x0) + SQR(y) + SQR(z) );
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
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

  // Go through cells
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {

        // Calculate position
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real r = get_r(x,y,z,x_exp);
        if (r >= r_m) {
          continue;
        }

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x - x_exp) / r * v;
	Real vy = y / r * v;
	Real vz = z / r * v;
	
	Real rho_val = get_rho_gauss(v,t_since_exp);

        if (rho_val <= rho_floor) {
          continue;
        }

        // Calculate energy
	// assume entropy ~ P/rho**(5/3) is constant
	// Get P from that
	Real u;
	if (GENERAL_EOS) {
	  Real p_val = Prho53*std::pow(rho_val,5./3.);
	  u = pmb->peos->EgasFromRhoP(rho_val, p_val);
	} else {
	  u = 3./2.*Prho53*std::pow(rho_val,5./3.);
	}
	    

        // Calculate effect of orbital velocity
        vy -= v_orb;
        Real ke = 0.5 * rho_val * (SQR(vx) + SQR(vy) + SQR(vz));

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);

	if (r < r_m and rho_val > rho_floor and vx > 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);
	
	if (r < r_m and rho_val > rho_floor and vx < 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);
	
	if (r < r_m and rho_val > rho_floor and vy > 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);
	
	if (r < r_m and rho_val > rho_floor and vy < 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);

	if (r < r_m and rho_val > rho_floor and vz > 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  AthenaArray<Real> &prim_scalar = pmb->pscalars->r;
  
  // Calculate ejecta parameters
  Real t_since_exp = time - t_exp;
  Real r_m = v_m * t_since_exp;
  
  Real rho_floor = pmb->peos->GetDensityFloor();
  // Real gamma = pmb->peos->GetGamma();

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
        Real r = get_r(x(i),y(j),z(k),x_exp);

	// Calculate velocity
	Real v = r / t_since_exp;
	Real vx = (x(i) - x_exp) / r * v;
	Real vy = y(j) / r * v;
	Real vz = z(k) / r * v;

	Real rho_val = get_rho_gauss(v,t_since_exp);

	if (r < r_m and rho_val > rho_floor and vz < 0.0) {
	  
	  // Calculate pressure
	  // assume entropy ~ P/rho**(5/3) is constant
	  // Get P from that
	  Real p_val = Prho53*std::pow(rho_val,5./3.);

	  // Calculate effect of orbital velocity
	  vy -= v_orb;

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
  return;
}

//----------------------------------------------------------------------------------------
// Refinement condition
// Inputs:
//   pmb: pointer to MeshBlock
// Outputs:
//   returned value: 1 (refine), -1 (derefine), or 0 (maintain level)
// Notes:
//   Checks gradient in composition fraction.

int GradientRefine(MeshBlock *pmb) {

  // Prepare index bounds
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Extract composition
  AthenaArray<Real> comp;
  comp.InitWithShallowSlice(pmb->pscalars->r, 4, 0, 1);

  // Count cells exceeding limits
  int num_cells = (ke - ks + 1) * (je - js + 1) * (ie - is + 1);
  int num_cells_refine = 0;
  int num_cells_derefine = 0;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real grad_x = (comp(0,k,j,i+1) - comp(0,k,j,i-1)) / 2.0;
        Real grad_y = (comp(0,k,j+1,i) - comp(0,k,j-1,i)) / 2.0;
        Real grad_z = (comp(0,k+1,j,i) - comp(0,k-1,j,i)) / 2.0;
        Real grad = grad_x;
        grad = std::max(grad, grad_y);
        grad = std::max(grad, grad_z);
        if (grad > refine_grad_val) {
          num_cells_refine++;
        }
        if (grad < derefine_grad_val) {
          num_cells_derefine++;
        }
      }
    }
  }

  // Evaluate refinement
  Real f_refine = static_cast<Real>(num_cells_refine) / static_cast<Real>(num_cells);
  Real f_derefine = static_cast<Real>(num_cells_derefine) / static_cast<Real>(num_cells);
  if (f_refine > refine_grad_frac) {
    return 1;
  }
  if (f_derefine > derefine_grad_frac) {
    return -1;
  }
  return 0;
}

// refinement condition: check the maximum pressure gradient
// and maximum scalar gradient.
// refine only if donor scalar is greater than r0_refine. 
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> &r = pmb->pscalars->r;
  Real maxeps = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
	if (r(0,k,j,i) >= r0_refine) {
	  Real epsp = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
				+SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
				+SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
	  Real epsr = std::sqrt(SQR(0.5*(r(0,k,j,i+1) - r(0,k,j,i-1)))
				+ SQR(0.5*(r(0,k,j+1,i) - r(0,k,j-1,i)))
				+ SQR(0.5*(r(0,k+1,j,i) - r(0,k-1,j,i))));
	  Real eps = std::max(epsp,epsr);
	  maxeps = std::max(maxeps, eps);
	}
      }
    }
  }

  if (maxeps > maxeps_threshold) return 1;
  if (maxeps < 0.25*maxeps_threshold) return -1;
  return 0;
}

//----------------------------------------------------------------------------------------
// Function called once every time step for user-defined work.
// Inputs: (none)
// Outputs: (none)
//void Mesh::UserWorkInLoop() {
  //MeshBlock *pmb = my_blocks(0);
  //Real bulk_mass = 0.0;
  //Real vx_integral = 0.0, vy_integral = 0.0, vz_integral = 0.0;
  //Real x_integral = 0.0, y_integral = 0.0, z_integral = 0.0;

  //for (int b = 0; b < nblocal; ++b) {
    //pmb = my_blocks(b);
    //Hydro *phydro = pmb->phydro;
    //Coordinates *pcoord = pmb->pcoord;
    //Gravity *pgrav = pmb->pgrav;
    //PassiveScalars *psclr = pmb->pscalars;

    // extract indices
    //int is = pmb->is, ie = pmb->ie;
    //int js = pmb->js, je = pmb->je;
    //int ks = pmb->ks, ke = pmb->ke;

    //Real enthalpy;

    // loop over element
    //for (int k = ks; k <= ke; ++ k){
      //for (int j = js; j <= je; ++ j){
	//for (int i = is; i <= ie; ++ i){
	  //Real rho = phydro->u(IDN,k,j,i);
	  //Real dV = pcoord->GetCellVolume(k,j,i);
	  //Real phi = pgrav->phi(k,j,i);
	  //Real s = psclr->s(0,k,j,i);

	  //if (time <= t_ej) {
	    // define star as where donor concentration is greater than 0.5
	    // if cell is bound and slow, add to mass
	    //if (phi <= phi_thres && s/rho >= r0_thres) {
	      //bulk_mass += rho * dV;
	      //vx_integral += phydro->u(IM1,k,j,i) * dV;
	      //vy_integral += phydro->u(IM2,k,j,i) * dV;
	      //vz_integral += phydro->u(IM3,k,j,i) * dV;
	      //x_integral += pcoord->x1v(i) * rho * dV;
	      //y_integral += pcoord->x2v(j) * rho * dV;
	      //z_integral += pcoord->x3v(k) * rho * dV;
	    //}
	  //} else {
	    // use Bernoulli parameter to define bound
	    // use bulk v from previous timestep to define 
	    //Real IMx = phydro->u(IM1,k,j,i);
	    //Real IMy = phydro->u(IM2,k,j,i);
	    //Real IMz = phydro->u(IM3,k,j,i);
	    //Real press = phydro->w(IPR,k,j,i);
	    
	    //Real ke = 0.5 * (SQR(IMx/rho - bulk_vx) + SQR(IMy/rho - bulk_vy) + SQR(IMz/rho - bulk_vz));

	    //if (GENERAL_EOS) {
	      //enthalpy = (pmb->peos->EgasFromRhoP(rho, press) + press)/rho;
	    //} else {
	      //Real gamma = pmb->peos->GetGamma();
	      //enthalpy = gamma/(gamma - 1.0) * press/rho;
	    //}
	    
	    //Real Be = ke + enthalpy + phi;

	    //if (Be < 0.0) {
	      //bulk_mass += rho * dV;
	      //vx_integral += phydro->u(IM1,k,j,i) * dV;
	      //vy_integral += phydro->u(IM2,k,j,i) * dV;
	      //vz_integral += phydro->u(IM3,k,j,i) * dV;
	      //x_integral += pcoord->x1v(i) * rho * dV;
	      //y_integral += pcoord->x2v(j) * rho * dV;
	      //z_integral += pcoord->x3v(k) * rho * dV;
	    //}
	  //}
	//}
      //}
    //}
  //}
//#ifdef MPI_PARALLEL
  //{
  //MPI_Allreduce(MPI_IN_PLACE, &bulk_mass, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &vx_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &vy_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &vz_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &x_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &y_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &z_integral, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  //}
//#endif

  // get mass-averaged velocity
  //bulk_vx = vx_integral / bulk_mass;
  //bulk_vy = vy_integral / bulk_mass;
  //bulk_vz = vz_integral / bulk_mass;

  // get C.O.M. of donor
  //x_star = x_integral / bulk_mass;
  //y_star = y_integral / bulk_mass;
  //z_star = z_integral / bulk_mass;

  // std::cout << bulk_mass << std::endl;
  // std::cout << "x_integral is" << x_integral << std::endl;
  // if (Globals::my_rank == 0) {
  //   std::cout << bulk_mass << std::endl;
  //   MPI_Bcast(&bulk_vx, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  //   MPI_Bcast(&bulk_vy, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  //   MPI_Bcast(&bulk_vz, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  //   std::cout << bulk_vx << std::endl;
  //   std::cout << bulk_vy << std::endl;
  //   std::cout << bulk_vz << std::endl;
  //   std::cout << x_star << std::endl;
  //   std::cout << y_star << std::endl;
  //   std::cout << z_star << std::endl;
  // }

  // if (time >= t_freeze) {
    
  // }
  
  //return;
//}

//----------------------------------------------------------------------------------------
// Get center of mass of donor defined as Phi < phi_thres and r0 > r0_thres

Real GetDonorCOM(MeshBlock *pmb, int iout) {
  if (iout == 0) {
    // std::cout << x_star << std::endl;
    return x_star;
  } else if (iout == 1) {
    return y_star;
  } else if (iout == 2) {
    return z_star;
  } else {
    return 0.0;
  }
}


//----------------------------------------------------------------------------------------
// Get momentum of donor defined as Phi < phi_thres and r0 > r0_thres
Real GetBulkV(MeshBlock *pmb, int iout) {
  if (iout == 3) {
    return bulk_vx;
  } else if (iout == 4) {
    return bulk_vy;
  } else if (iout == 5) {
    return bulk_vz;
  } else {
    return 0.0;
  }
}

//----------------------------------------------------------------------------------------
// Get total bound mass, bound mass of r0 material, and bound mass of r1 material
// bound defined as Be < 0

//Comment out Bound mass - Kathlynn? 
Real GetBoundMass(MeshBlock *pmb, int iout) {
  Real bound_mass = 0.0;

  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  Gravity *pgrav = pmb->pgrav;
  PassiveScalars *psclr = pmb->pscalars;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	Real IMx = phydro->u(IM1,k,j,i);
	Real IMy = phydro->u(IM2,k,j,i);
	Real IMz = phydro->u(IM3,k,j,i);
	Real rho = phydro->u(IDN,k,j,i);
	Real press = phydro->w(IPR,k,j,i);
	Real dV = pcoord->GetCellVolume(k,j,i);
	Real phi = pgrav->phi(k,j,i);
	Real s = 0.0;

	// calculate Bernoulli parameter relative to bulk motion
	Real ke = 0.5 * (SQR(IMx/rho - bulk_vx) + SQR(IMy/rho - bulk_vy) + SQR(IMz/rho - bulk_vz));

	Real enthalpy;
	if (GENERAL_EOS) {
	  enthalpy = (pmb->peos->EgasFromRhoP(rho, press) + press)/rho;
	} else {
	  Real gamma  = pmb->peos->GetGamma();
	  enthalpy = gamma/(gamma - 1.0) * press/rho;
	}
	
	Real Be = ke + enthalpy + phi;

	// if Bernoulli parameter is negative, treat as bound
	if (Be < 0.0) {
	  if (iout == 6) {
	    bound_mass += rho * dV;
	  } else if (iout == 7) {
	    s = psclr->s(0,k,j,i);
	    bound_mass += s * dV;
	  } else if (iout == 8) {
	    s = psclr->s(1,k,j,i);
	    bound_mass += s * dV;
	  }
	}
      }
    }
  }
  return bound_mass;
}

//----------------------------------------------------------------------------------------
// Get momentum of bound material
// bound defined as Be < 0

Real GetBoundMom(MeshBlock *pmb, int iout) {
  Real bound_v = 0.0;

  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  Gravity *pgrav = pmb->pgrav;

  //std::cout << bulk_vx << std::endl;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	Real IMx = phydro->u(IM1,k,j,i);
	Real IMy = phydro->u(IM2,k,j,i);
	Real IMz = phydro->u(IM3,k,j,i);
	Real rho = phydro->u(IDN,k,j,i);
	Real press = phydro->w(IPR,k,j,i);
	Real dV = pcoord->GetCellVolume(k,j,i);
	Real phi = pgrav->phi(k,j,i);

	// calculate Bernoulli parameter relative to bulk motion
	Real ke = 0.5 * (SQR(IMx/rho - bulk_vx) + SQR(IMy/rho - bulk_vy) + SQR(IMz/rho - bulk_vz));
	
	Real enthalpy;
	if (GENERAL_EOS) {
	  enthalpy = (pmb->peos->EgasFromRhoP(rho, press) + press)/rho;
	} else {
	  Real gamma  = pmb->peos->GetGamma();
	  enthalpy = gamma/(gamma - 1.0) * press/rho;
	}
	
	Real Be = ke + enthalpy + phi;

	// if Bernoulli parameter is negative, treat as bound
	if (Be < 0.0) {
	  if (iout == 9) {
	    bound_v += IMx * dV;
	  } else if (iout == 10) {
	    bound_v += IMy * dV;
	  } else if (iout == 11) {
	    bound_v += IMz * dV;
	  }
	}
      }
    }
  }
  return bound_v;
}

//----------------------------------------------------------------------------------------
// Get central pressure/density at the donor C.O.M.

Real GetDonorPDc(MeshBlock *pmb, int iout) {
  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  Gravity *pgrav = pmb->pgrav;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real PD_c = 0.0;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	// get cell boundaries
	Real x1l = pcoord->x1f(i), x1u = pcoord->x1f(i+1);
	Real x2l = pcoord->x2f(j), x2u = pcoord->x2f(j+1);
	Real x3l = pcoord->x3f(k), x3u = pcoord->x3f(k+1);

	// if C.O.M. is inside of cell, get quantity
	if (x1l <= x_star && x_star <= x1u && x2l <= y_star && y_star <= x2u && x3l <= z_star && z_star <= x3u) {
	  if (iout == 12) {
	    PD_c = phydro->u(IDN,k,j,i);
	  } else if (iout == 13) {
	    PD_c = phydro->w(IPR,k,j,i);
	  }
	}
      }
    }
  }
  return PD_c;
}

//----------------------------------------------------------------------------------------
// Get max density/pressure inside donor

Real GetMaxPD(MeshBlock *pmb, int iout) {
  Hydro *phydro = pmb->phydro;
  PassiveScalars *psclr = pmb->pscalars;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real max_PD = 0.0;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	Real s = psclr->s(0,k,j,i);
	Real den = phydro->u(IDN,k,j,i);
	// get quantity only if donor concentration is high enough
	if (s/den > r0_thres) {
	  if (iout == 14) {
	    max_PD = std::max(max_PD,phydro->u(IDN,k,j,i));
	  } else if (iout == 15) {
	    max_PD = std::max(max_PD,phydro->w(IPR,k,j,i));
	  }
	}
      }
    }
  }
  return max_PD;
}

//----------------------------------------------------------------------------------------
// get total mass at boundary

Real GetMdotBound(MeshBlock *pmb, int iout) {
  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  PassiveScalars *psclr = pmb->pscalars;
  Gravity *pgrav = pmb->pgrav;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real myval = 0.0;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	// get cell boundaries
	Real x1l = pcoord->x1f(i), x1u = pcoord->x1f(i+1);
	Real x2l = pcoord->x2f(j), x2u = pcoord->x2f(j+1);
	Real x3l = pcoord->x3f(k), x3u = pcoord->x3f(k+1);

	Real rho = phydro->u(IDN,k,j,i);

	Real vx = phydro->u(IM1,k,j,i) / rho;
	Real vy = phydro->u(IM2,k,j,i) / rho;
	Real vz = phydro->u(IM3,k,j,i) / rho;

	Real s = psclr->s(0,k,j,i);
	  
	Real press = phydro->w(IPR,k,j,i);
	Real phi = pgrav->phi(k,j,i);

	// calculate Bernoulli parameter relative to bulk motion
	Real ke = 0.5 * (SQR(vx - bulk_vx) + SQR(vy - bulk_vy) + SQR(vz - bulk_vz));
	  
	Real enthalpy;
	if (GENERAL_EOS) {
	  enthalpy = (pmb->peos->EgasFromRhoP(rho, press) + press)/rho;
	} else {
	  Real gamma  = pmb->peos->GetGamma();
	  enthalpy = gamma/(gamma - 1.0) * press/rho;
	}
	
	Real Be = ke + enthalpy + phi;
	  
	if (Be <= 0.0) {
	  
	  // left X face
	  if (x1l == pmb->pmy_mesh->mesh_size.x1min and vx < 0.0) {
	    Real dA = pcoord->GetFace1Area(k,j,i);
	    myval += s * std::abs(vx) * dA;
	    // right X face
	  } else if (x1u == pmb->pmy_mesh->mesh_size.x1max and vx > 0.0) {
	    Real dA = pcoord->GetFace1Area(k,j,i);
	    myval += s * std::abs(vx) * dA;
	    // front Y face 
	  } else if (x2l == pmb->pmy_mesh->mesh_size.x2min and vy < 0.0) {
	    Real dA = pcoord->GetFace2Area(k,j,i);
	    myval += s * std::abs(vy) * dA;
	    // back Y face 
	  } else if (x2u == pmb->pmy_mesh->mesh_size.x2max and vy > 0.0) {
	    Real dA = pcoord->GetFace2Area(k,j,i);
	    myval += s * std::abs(vy) * dA;
	    // bot Z face 
	  } else if (x3l == pmb->pmy_mesh->mesh_size.x3min and vz < 0.0) {
	    Real dA = pcoord->GetFace3Area(k,j,i);
	    myval += s * std::abs(vz) * dA;
	    // top Z face 
	  } else if (x3u == pmb->pmy_mesh->mesh_size.x3max and vz > 0.0) {
	    Real dA = pcoord->GetFace3Area(k,j,i);
	    myval += s * std::abs(vz) * dA;
	  }
	}
      }
    }
  }
	
  return myval;
}


//----------------------------------------------------------------------------------------
// get total mass at boundary
Real GetM_v(MeshBlock *pmb, int iout) {
  Hydro *phydro = pmb->phydro;
  Coordinates *pcoord = pmb->pcoord;
  PassiveScalars *psclr = pmb->pscalars;
  Gravity *pgrav = pmb->pgrav;

  // extract indices
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real myval = 0.0;

  // loop over element
  for (int k = ks; k <= ke; ++ k){
    for (int j = js; j <= je; ++ j){
      for (int i = is; i <= ie; ++ i){
	// get cell boundaries
	Real x1l = pcoord->x1f(i), x1u = pcoord->x1f(i+1);
	Real x2l = pcoord->x2f(j), x2u = pcoord->x2f(j+1);
	Real x3l = pcoord->x3f(k), x3u = pcoord->x3f(k+1);

	// get angle from explosion center
	Real x = pcoord->x1v(i) - x_exp;
	Real y = pcoord->x2v(j);
	Real z = pcoord->x3v(k);
	Real r_yz = std::sqrt( y*y + z*z );
	
	Real theta = std::atan(std::abs(r_yz/x)) * 180.0/PI;
	
	// get other quantities

	Real rho = phydro->u(IDN,k,j,i);

	Real vx = phydro->u(IM1,k,j,i) / rho;
	Real vy = phydro->u(IM2,k,j,i) / rho;
	Real vz = phydro->u(IM3,k,j,i) / rho;

	Real vmag = std::sqrt( SQR(vx) + SQR(vy) + SQR(vz) )*v_0;
	Real logvmag = std::log10( vmag );
	bool in_bin = false;

	Real s = 0.0;
	int sclr_index = 0; // donor material by default

	// assign material and check if cell is in v bin or cos(theta) bin
	if (iout < 63) {
	  // donor material
	  if (iout == 16 and logvmag > 9.3) {
	    in_bin = true;
	  } else {
	    Real logvmin = 9.3 - 0.05*(iout-16);
	    Real logvmax = logvmin + 0.05;
	    if (logvmag <= logvmax && logvmag > logvmin) {
	      in_bin = true;
	    }
	  }
	} else if (iout >= 63 && iout < 111) {
	  // SN material 
	  sclr_index = 1;
	  
	  if (iout == 63 and logvmag > 9.3) {
	    in_bin = true;
	  } else {
	    Real logvmin = 9.3 - 0.05*(iout-63);
	    Real logvmax = logvmin + 0.05;
	    if (logvmag <= logvmax && logvmag > logvmin) {
	      in_bin = true;
	    }
	  }
	} else if (iout >= 111 && iout < 156) {
	  // donor material

	  Real theta_min = 2.0*(iout-111);
	  Real theta_max = theta_min + 2.0;

	  if (theta < theta_max && theta >= theta_min) {
	    in_bin = true;
	  }
	} else {
	  // SN ejecta material
	  sclr_index = 1;

	  Real theta_min = 2.0*(iout-156);
	  Real theta_max = theta_min + 2.0;

	  if (theta < theta_max && theta >= theta_min) {
	    in_bin = true;
	  }
	}

	s = psclr->s(sclr_index,k,j,i);
	
	if (in_bin) {
	  
	  Real press = phydro->w(IPR,k,j,i);
	  Real phi = pgrav->phi(k,j,i);

	  // calculate Bernoulli parameter relative to bulk motion
	  Real ke = 0.5 * (SQR(vx - bulk_vx) + SQR(vy - bulk_vy) + SQR(vz - bulk_vz));

	  Real enthalpy;
	  if (GENERAL_EOS) {
	    enthalpy = (pmb->peos->EgasFromRhoP(rho, press) + press)/rho;
	  } else {
	    Real gamma  = pmb->peos->GetGamma();
	    enthalpy = gamma/(gamma - 1.0) * press/rho;
	  }
	
	  Real Be = ke + enthalpy + phi;

	  if (Be > 0.0) {
	  
	    // left X face
	    if (x1l == pmb->pmy_mesh->mesh_size.x1min and vx < 0.0) {
	      Real dA = pcoord->GetFace1Area(k,j,i);
	      // do not count SN material in case of fall back
	      if (sclr_index == 0) {
		myval += s * std::abs(vx) * dA;
	      }
	      // right X face
	    } else if (x1u == pmb->pmy_mesh->mesh_size.x1max and vx > 0.0) {
	      Real dA = pcoord->GetFace1Area(k,j,i);
	      myval += s * std::abs(vx) * dA;
	      // front Y face 
	    } else if (x2l == pmb->pmy_mesh->mesh_size.x2min and vy < 0.0) {
	      Real dA = pcoord->GetFace2Area(k,j,i);
	      myval += s * std::abs(vy) * dA;
	      // back Y face 
	    } else if (x2u == pmb->pmy_mesh->mesh_size.x2max and vy > 0.0) {
	      Real dA = pcoord->GetFace2Area(k,j,i);
	      myval += s * std::abs(vy) * dA;
	      // bot Z face 
	    } else if (x3l == pmb->pmy_mesh->mesh_size.x3min and vz < 0.0) {
	      Real dA = pcoord->GetFace3Area(k,j,i);
	      myval += s * std::abs(vz) * dA;
	      // top Z face 
	    } else if (x3u == pmb->pmy_mesh->mesh_size.x3max and vz > 0.0) {
	      Real dA = pcoord->GetFace3Area(k,j,i);
	      myval += s * std::abs(vz) * dA;
	    }
	  }
	}
      }
    }
  }
	
  return myval;
}

