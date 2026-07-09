//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cee_sn.cpp
//! \brief Initializes CEE outflows within domain and sets SN inner radial boundary.

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

//#include <unistd.h>   
#include <stdio.h>

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

// inflow/outflow BCs
void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// vacuum boundary
void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh);

Real MyTimeStep(MeshBlock* pmb);

void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt, 
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar);

namespace {
std::vector<Real> vr_in, rho_in, temp_in, ar36_in, fe56_in, co56_in, ni56_in;
static Real vr_in_current, rho_in_current, temp_in_current, 
	    ar36_in_current, fe56_in_current, co56_in_current, ni56_in_current;
Real gammagas, Rgas, vmax, ramPressureFactor, rhoISM, r_inner, Rsun, Mej, Eej, t0;
Real Mdotwind, vwind;
Real day;
Real epsilon_Ni, epsilon_Co;
Real tau_Ni, tau_Co;
Real A_nuc, mproton;
Real mu_SN_ejecta, ar, kB;
bool diode;
int NumToRead;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gammagas          = pin->GetOrAddReal("hydro","gamma",0.0);
  diode             = pin->GetOrAddBoolean("problem","diode",false);
  vmax              = pin->GetOrAddReal("problem","v_max",0.0);
  ramPressureFactor = pin->GetOrAddReal("problem","ramPressureFactor",0.0);
  rhoISM            = pin->GetOrAddReal("problem","rho_ISM",0.0);
  r_inner           = pin->GetOrAddReal("mesh","x1min",0.0);
  Mej               = pin->GetOrAddReal("problem","Mej",0.0);
  Eej               = pin->GetOrAddReal("problem","Eej",0.0);
  Rsun              = 7.0e10;
  Rgas              = 8.314e7;
  Mdotwind          = 1.0e-8*2.0e33/365.25/24.0/3600.0;
  vwind             = 30.0e5;
  t0                = r_inner/vmax;
  day               = 24.0*3600.0;
  mproton           = 1.6726e-24;
  mu_SN_ejecta      = 2.0;
  ar                = 7.5646e-15; // radiation density constant
  kB                = 1.3807e-16; // Boltzmann constant
  NumToRead         = 107;
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiodeOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);
  EnrollUserTimeStepFunction(MyTimeStep);
  EnrollUserExplicitSourceFunction(radioactiveHeating);

  epsilon_Ni = 1.72e-6; // energy released in decays
  epsilon_Co = 3.49e-6;
  tau_Ni     = 8.77*day; // e-folding times of decays
  tau_Co     = 111.0*day;
  A_nuc      = 56.0; // mass number
  
  if (Globals::my_rank==0) {
  char vrFile[256], rhoFile[256], tempFile[256], 
       ar36File[256], fe56File[256], co56File[256], ni56File[256];
  sprintf(vrFile,   "athenainput_vr.txt");
  sprintf(rhoFile,  "athenainput_rho.txt");
  sprintf(tempFile, "athenainput_temp.txt");
  sprintf(ar36File, "athenainput_ar36.txt");
  sprintf(fe56File, "athenainput_fe56.txt");
  sprintf(co56File, "athenainput_co56.txt");
  sprintf(ni56File, "athenainput_ni56.txt");
  printf("Opening data files with state variables...\n");
  std::ifstream vrFileRead, rhoFileRead, tempFileRead, 
                ar36FileRead, fe56FileRead, co56FileRead, ni56FileRead;
  vrFileRead.open(vrFile);
  rhoFileRead.open(rhoFile);
  tempFileRead.open(tempFile);
  ar36FileRead.open(ar36File);
  fe56FileRead.open(fe56File);
  co56FileRead.open(co56File);
  ni56FileRead.open(ni56File);

  Real vr, rho, temp, ar36, fe56, co56, ni56;
  for (int l = 0; l < NumToRead; l++) {
    vrFileRead   >> vr;
    rhoFileRead  >> rho;
    tempFileRead >> temp;
    ar36FileRead >> ar36;
    fe56FileRead >> fe56;
    co56FileRead >> co56;
    ni56FileRead >> ni56;
    vr_in.push_back(vr);
    rho_in.push_back(rho);
    temp_in.push_back(temp);
    ar36_in.push_back(ar36);
    fe56_in.push_back(fe56);
    co56_in.push_back(co56);
    ni56_in.push_back(ni56);
  }
  printf("Done reading, closing data files\n");
  vrFileRead.close();
  rhoFileRead.close();
  tempFileRead.close();
  ar36FileRead.close();
  fe56FileRead.close();
  co56FileRead.close();
  ni56FileRead.close();

  rho_in_current  = rho_in[NumToRead-1];
  vr_in_current   = vr_in[NumToRead-1];
  temp_in_current = temp_in[NumToRead-1];
  ar36_in_current = ar36_in[NumToRead-1];
  fe56_in_current = fe56_in[NumToRead-1];
  co56_in_current = co56_in[NumToRead-1];
  ni56_in_current = ni56_in[NumToRead-1];
  }
  
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes CEE outflows within domain.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real r, theta, z;
  Real diskHeight, rhoCEE, rhoMin, rhoWind;
  Real rho, temp, pres, mintemp;

  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    //phi = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      theta = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);
        z = r*std::cos(theta);

        diskHeight = 0.4*r+10.0*Rsun;
		// Rsun*(95.0*std::log10(r/Rsun)-125.0);
	rhoCEE = std::exp(-z*z/2.0/diskHeight/diskHeight)
	//* 0.01*std::pow(r/10.0/Rsun,-4.0)*std::pow(1.0+std::pow(125.0*Rsun/r,3.5),-1.05); // 500 d
	* (1.0e-4*std::pow(r/100.0/Rsun,-4.8)*std::pow(1.0+std::pow(690.0*Rsun/r,1.9),-3.4)+6.358e-13*std::pow(r/6000.0/Rsun,-3.2)/(1.0+std::pow(6000.0*Rsun/r,6.0))); // 10000 d
	
	//rhoMin = 1.0e-11;
	rhoWind = Mdotwind/4.0/3.14159/r/r/vwind;

	mintemp = 1.0e4;

	if (rhoCEE > rhoWind && rhoCEE > rhoISM) { // disk
	  rho = rhoCEE;
    temp = std::max( mintemp, 4.5e4/(r/100.0/Rsun) );
    pres = rho*Rgas*temp;
	} else if (rhoWind > rhoISM) { // wind
	  rho = rhoWind;
    pres = rho*Rgas*mintemp;
	} else { // ISM
    rho = rhoISM;
    pres = ramPressureFactor*rhoISM*vmax*vmax;
  }

        phydro->u(IDN,k,j,i) = rho;

        phydro->u(IM1,k,j,i) = 0.0; // rho*vel0*std::cos(x2); // radial
        phydro->u(IM2,k,j,i) = 0.0; //-rho*vel0*std::sin(x2); // polar
        phydro->u(IM3,k,j,i) = 0.0;               // azimuth

        phydro->u(IEN,k,j,i) = pres;
		            // ramPressureFactor*rhoCEE*vmax*vmax; 
                            // pres/(gammagas-1.0) + 0.5*rho*vel0*vel0;

        for (int l=0; l<NSCALARS; ++l) {
          pscalars->s(l,k,j,i) = 0.0;
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void DiodeOuterX1()
//  \brief Sets doide outflow conditions at outer x1 boundary
//
// Gas outflows with optional diode condition

void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
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
	
        for (int l=0; l<NSCALARS; ++l) {
          pmb->pscalars->r(l,k,j,iu+i) = pmb->pscalars->r(l,k,j,iu);
        }

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
//  \brief Sets supernova ejecta inner boundary
//
// Quantities in ghost cells are set to Gaussian ejecta model from Wong+24

void SNInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt,
               int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

  Real v_inner, pres;
  Real rhoSunny, v0sq, prefactor;

  v_inner = r_inner/(time + t0);

  v0sq = 4.0 / 3.0 * Eej / Mej;
  prefactor = std::pow(3.0 / 4.0 / 3.14159 / Eej, 1.5) * std::pow(Mej, 2.5);
  rhoSunny = prefactor * std::exp(-v_inner * v_inner / v0sq) * std::pow(time + t0, -3.0);
  pres = 0.7e14*std::pow(rhoSunny,1.6666666666667);
  //pres = ramPressureFactor*rhoSunny*vmax*vmax;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {

      	//pres = rho_in_current*Rgas*temp_in_current/mu_SN_ejecta
        //     + ar*std::pow(temp_in_current,3)/3.0;

        prim(IDN,k,j,il-i)        = rhoSunny; // rho_in_current;
        prim(IM1,k,j,il-i)        = v_inner;
        prim(IM2,k,j,il-i)        = 0.0;
        prim(IM3,k,j,il-i)        = 0.0;
        prim(IEN,k,j,il-i)        = pres;
	/*
        pmb->pscalars->r(0,k,j,il-i) = 1.0;
        pmb->pscalars->r(1,k,j,il-i) = 0.0; // ar36_in_current;
        pmb->pscalars->r(2,k,j,il-i) = 0.0; // fe56_in_current;
        pmb->pscalars->r(3,k,j,il-i) = 0.0; // co56_in_current;
        pmb->pscalars->r(4,k,j,il-i) = 0.0; // ni56_in_current;
	*/
      }
    }
  }
}

Real MyTimeStep(MeshBlock* pmb) {
    Real time   = pmb->pmy_mesh->time;
    Real dt     = pmb->pmy_mesh->dt;
    Real min_dt = 1.0e3;
    if (time < 0.1) min_dt = std::min(dt, 1.0e-2);
    return min_dt;
}

void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt,
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar) {

    Real deltaRho_Ni, deltaRho_Co;

    Real coefficient_Ni = epsilon_Ni / (tau_Ni * A_nuc * mproton);
    Real coefficient_Co = epsilon_Co / (tau_Co * A_nuc * mproton);

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
            for (int i = pmb->is; i <= pmb->ie; ++i) {
                cons(IEN,k,j,i) += dt * cons_scalar(4,k,j,i) * coefficient_Ni 
                                 * std::exp(-1.0 * (time + t0) / tau_Ni)
                                 + dt * cons_scalar(3,k,j,i) * coefficient_Co
                                 * std::exp(-1.0 * (time + t0) / tau_Co);
                deltaRho_Ni = -dt/tau_Ni*cons_scalar(4,k,j,i); // Ni to Co
                deltaRho_Co = -dt/tau_Co*cons_scalar(3,k,j,i); // Co to Fe
                cons_scalar(4,k,j,i) += deltaRho_Ni;
                cons_scalar(3,k,j,i) += deltaRho_Co - deltaRho_Ni;
                cons_scalar(2,k,j,i) += -deltaRho_Co;
            }
        }
    }
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void Mesh::UserWorkInLoop() {

  if (Globals::my_rank==0) {

    Real v_inner = r_inner/(time + t0);

    Real dist;
    int index;

    dist = 1.0e10;
    index = -1;

    //std::cout << "starting search" << std::endl;
    for (int l=0; l<NumToRead; ++l) {
      if (dist > std::abs(v_inner-vr_in[l])) {
        dist = std::abs(v_inner-vr_in[l]);
        index = l;
      } else {
        break;
      }
    }
    std::cout << "vr inner = " << v_inner << std::endl;
    std::cout << "Found index " << index << " with vr = " << vr_in[index] << " and rho = " << rho_in[index] << std::endl;

    rho_in_current  = rho_in[index];
    vr_in_current   = vr_in[index];
    temp_in_current = temp_in[index];
    ar36_in_current = ar36_in[index];
    fe56_in_current = fe56_in[index];
    co56_in_current = co56_in[index];
    ni56_in_current = ni56_in[index];

  }

}

namespace {
}
