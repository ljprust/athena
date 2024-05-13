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
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cfloat>

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
void SNOuterX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// inject SN ejecta
void SNInnerX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
Real MyTimeStep(MeshBlock* pmb);
namespace {
    // problem parameters which are useful to make global to this file
    Real rho0, p0, vel0, gammagas, pressEj,rhoEj, minRadius;
    bool diode;
    Real Rej, Eej, Mej, vmax;
} // namespace
//Some Useful Functions
//gets pressure given temperature and density
Real pressureEqPgen(Real temp, Real rho) {

  Real a=7.56*std::pow(10.0,-15.0);
  Real mu=2.0;
  Real mProton=1.6726*std::pow(10.0,-24.0);
  Real kB=1.3807*std::pow(10.0,-16.0);
  return ((1/3.)*a*std::pow(temp,4.0))+rho*kB*temp/(mu*mProton);
}

//gets internal energy given temperature and density
Real energyEqPgen(Real temp, Real rho) {

  Real a=7.56*std::pow(10.0,-15.0);
  Real mu=2.0;
  Real mProton=1.6726*std::pow(10.0,-24.0);
  Real kB=1.3807*std::pow(10.0,-16.0);
  return (a*std::pow(temp,4.0))+1.5* rho * kB * temp / (mu * mProton);
}

Real findRoot(Real A, Real B) {
    Real z3 = std::pow(81.0 * std::pow(A, 4.0) + 768.0 * std::pow(B, 3.0), .5) + 9.0 * A * A;
    Real numerator = std::pow(2.0, 1 / 3.) * std::pow(z3, 2.0 / 3.0) - 8.0 * std::pow(3.0, 1.0 / 3.0) * B;
    Real denominator = std::pow(6.0, 2.0 / 3.0) * std::pow(z3, 1.0 / 3.0);
    return numerator / denominator;

}
//takes in pressure and density and solves the quartic analytically to get you temperature
Real calcTemperaturePressure(Real rho, Real pres) {
    Real a, mu, mProton, kB, A, B, temp, z1, z2, y;
    mProton = 1.6726 * std::pow(10.0, -24.0);
    kB = 1.3807 * std::pow(10.0, -16.0);
    mu = 2.0;
    a = 7.56 * std::pow(10.0, -15.0);
    A = 3.0 * kB * rho / (a * mu * mProton);
    B = 3.0 * pres / a;
    y = findRoot(A, B);
    temp = std::pow(y, .5) * (std::pow(2.0 * A / std::pow(y * y * y, .5) - 1, .5) - 1) / 2;
    return temp;
}
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

// read parameters from input file
void Mesh::InitUserMeshData(ParameterInput* pin) {
    // ambient medium properties
    rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
    p0 = pin->GetOrAddReal("problem", "p0", 1.0);
    vel0= pin->GetOrAddReal("problem","vel0",1.0);
    // gamma law for gas
    gammagas = pin->GetOrAddReal("hydro", "gamma", 1.66666667);
    // use diode condition for outflow
    diode = pin->GetOrAddBoolean("problem", "diode", false);
    // properties of ejecta
    Rej = pin->GetOrAddReal("problem", "Rej", 1.0);
    Eej = pin->GetOrAddReal("problem", "Eej", 1.0);
    Mej = pin->GetOrAddReal("problem", "Mej", 1.0);
    rhoEj=pin->GetOrAddReal("problem","rhoEj",1.0);
    pressEj=pin->GetOrAddReal("problem","pressEj",1.0);
    vmax = pin->GetOrAddReal("problem", "vmax", 1.0);
    minRadius = pin->GetOrAddReal("mesh", "x1min", 1.0);
        // tell athena to use our user-defined boundaries
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1);
     //tell athena to use user-defined timestep
    EnrollUserTimeStepFunction(MyTimeStep);
    return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes SN ejecta.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput* pin) {

    Real r, t0, v0sq, v, rho, pres, tChar, charEnergyDensity, temp;
    t0 = Rej / vmax;
    v0sq = 4.0 / 3.0 * Eej / Mej;
    //  Initialize density and momenta
    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {
                r = pcoord->x1v(i);

                // transverse velocities are always zero
                phydro->u(IM2, k, j, i) = 0.0;
                phydro->u(IM3, k, j, i) = 0.0;

                if (r > Rej) {
                    // ambient medium
                    temp=calcTemperaturePressurePgen(rho, p0);
                    phydro->u(IDN, k, j, i) = rho0;
                    phydro->u(IM1, k, j, i) =rho0*vel0;
                    phydro->u(IEN, k, j, i) = energyEqPgen(temp, rho0)+.5*rho0*vel0*vel0;
		            pscalars->s(0,k,j,i)=rho0;
                }
                else {
                    // ejecta profile
                    v = vmax*r/Rej;
                    rho = std::pow(3.0 / 4.0 / 3.14159 / Eej, 1.5) * std::pow(Mej, 2.5)
                        * std::exp(-v * v / v0sq) / t0 / t0 / t0;
                    pres =0.69*std::pow(10.0,14.0)*std::pow(rho,1.66666666666667);
                    temp=2*std::pow(10.0,5.0);
                    phydro->u(IDN, k, j, i) = rho;
                    phydro->u(IM1, k, j, i) =rho*v;
                    phydro->u(IEN, k, j, i) = 1.0*(.5)*rho*v*v+energyEqPgen(temp, rho);
		            pscalars->s(0,k,j,i)=0.0;
                }

            }
        }
    }

    return;
}

//----------------------------------------------------------------------------------------
//! \fn void SNOuterX1()
//  \brief Sets boundary condition on downstream outflow boundary

void SNOuterX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

    bool applyDiode;

    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {

                // let gas outflow freely
                prim(IDN, k, j, iu + i) = prim(IDN, k, j, iu);
                prim(IM2, k, j, iu + i) = prim(IM2, k, j, iu);
                prim(IM3, k, j, iu + i) = prim(IM3, k, j, iu);
                prim(IEN, k, j, iu + i) = prim(IEN, k, j, iu);

                // ensure that no gas enters through the outflow boundary
                // by giving the radial velocity some TLC
                applyDiode = prim(IM1, k, j, iu) < 0.0;
                if (diode && applyDiode) {
                    prim(IM1, k, j, iu + i) = 0.0;
                }
                else {
                    prim(IM1, k, j, iu + i) = prim(IM1, k, j, iu);
                }

            }
        }
    }
}


//--------------------------------------------------------------------------------
//Calculates the temperature of a thing given our pressure
//----------------------------------------------------------------------------------------
//! \fn void SNInnerX1()
//  \brief Sets inner boundary where SN ejecta is injected

void SNInnerX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

    Real rho, vel, v0sq, pres, prefactor, t0, tChar, charEnergyDensity, temp;
    t0 = Rej / vmax;
    v0sq = 4.0 / 3.0 * Eej / Mej;
    vel =minRadius/(time+t0);
    prefactor = std::pow(3.0 / 4.0 / 3.14159 / Eej, 1.5) * std::pow(Mej, 2.5);
    rho = prefactor * std::exp(-vel * vel / v0sq) * std::pow(time + t0, -3.0);
    pres =0.69*std::pow(10.0,14.0)*std::pow(rho,1.6666666666667);
    temp=2*std::pow(10.0,5.0);
    if (Rej == 0.0 || Eej == 0.0 || Mej == 0.0 || vmax == 0.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in sn.cpp ProblemGenerator" << std::endl
            << "ejecta parameters not set" << std::endl;
        ATHENA_ERROR(msg);
    }
  
    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {
                prim(IDN, k, j, il - i) = rho;
                prim(IM1, k, j, il - i) =vel;
                prim(IM2, k, j, il - i) = 0.0;
                prim(IM3, k, j, il - i) = 0.0;
                prim(IEN, k, j, il - i) = energyEqPgen(temp, rho);
            }
        }
    }
}

Real MyTimeStep(MeshBlock *pmb) {
  Real time = pmb ->pmy_mesh -> time;
  Real dt=pmb->pmy_mesh ->dt;
  Real min_dt=1e2;
  if (time<.01) {
          min_dt=std::min(dt,1e-4);	
  } // calculate your own time step here

  return min_dt;
}
/*
void SNInnerX1Wind(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

    Real rho, vel, v0sq, pres, prefactor, t0;

    vel = std::pow(10,5);
    rho = std::pow(10,4);
    pres =std::pow(10, 12)*3;

    if (Rej == 0.0 || Eej == 0.0 || Mej == 0.0 || vmax == 0.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in sn.cpp ProblemGenerator" << std::endl
            << "ejecta parameters not set" << std::endl;
        ATHENA_ERROR(msg);
    }

    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {

                prim(IDN, k, j, il - i) = rho;
                prim(IM1, k, j, il - i) = vel;
                prim(IM2, k, j, il - i) = 0.0;
                prim(IM3, k, j, il - i) = 0.0;
                prim(IEN, k, j, il - i) = pres;

            }
        }
    }
}
*/
 

//namespace {}
