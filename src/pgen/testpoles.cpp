//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file octree_test.cpp
//! \brief Problem generator for reading in data using an octree search

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include<unistd.h>   
#include <stdio.h>
#include <iostream>

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
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"

int NumToRead, nx1, NumToReadEdge;
Real x1min, x1max, deltax, gammagas;
std::vector<Real> r_in, theta_in, x_in, y_in, rho_in, vr_in, pres_in, edgeR_in, edgeRho_in, edgeVr_in, edgePres_in;

// outflow condition with diode
void SNOuterX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SNOuterX1Rad(MeshBlock* pmb, Coordinates* pco, NRRadiation* pnrrad, const AthenaArray<Real>& prim, FaceField& b, AthenaArray<Real>& ir,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SNInnerX1Rad(MeshBlock* pmb, Coordinates* pco, NRRadiation* pnrrad, const AthenaArray<Real>& prim, FaceField& b, AthenaArray<Real>& ir,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// inject SN ejecta
void SNInnerX1Wind(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt,
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar);

void opacityCalc(MeshBlock* pmb, AthenaArray<Real>& prim);

Real MyTimeStep(MeshBlock* pmb);
namespace {
    // problem parameters which are useful to make global to this file
    Real rho0, p0, vel0, gammaGas, pressEj, rhoEj, vEj, minRadius, massFraction, T0, density0, l0, tempInitial, opacity, mu, initialTime;
    bool diode;
    Real Rej, Eej, Mej, vmax;
} // namespace
//Some Useful Functions
//gets pressure given temperature and density
Real pressureEqPgen(Real temp, Real rho) {

    Real a = 7.56 * std::pow(10.0, -15.0);
    Real mu = 4.0 / 3.0;
    Real mProton = 1.6726 * std::pow(10.0, -24.0);
    Real kB = 1.3807 * std::pow(10.0, -16.0);
    return ((1 / 3.) * a * std::pow(temp, 4.0)) + rho * kB * temp / (mu * mProton);
}

//gets internal energy given temperature and density
Real energyEqPgen(Real temp, Real rho) {

    Real a = 7.56 * std::pow(10.0, -15.0);
    Real mu = 4.0 / 3.0;
    Real mProton = 1.6726 * std::pow(10.0, -24.0);
    Real kB = 1.3807 * std::pow(10.0, -16.0);
    return (a * std::pow(temp, 4.0)) + 1.5 * rho * kB * temp / (mu * mProton);
}

Real gasEnergyEq(Real rho, Real temp) {
    Real mu = 4.0 / 3.0;
    Real mProton = 1.6726 * std::pow(10.0, -24.0);
    Real kB = 1.3807 * std::pow(10.0, -16.0);
    return 1.5 * rho * kB * temp / (mu * mProton);
}
//Helper functions to calculate temperature:
Real findCubeRoot(Real A, Real B) {
    Real z3 = std::pow(81.0 * std::pow(A, 4.0) + 768.0 * std::pow(B, 3.0), .5) + 9.0 * A * A;
    Real numerator = std::pow(2.0, 1 / 3.) * std::pow(z3, 2.0 / 3.0) - 8.0 * std::pow(3.0, 1.0 / 3.0) * B;
    Real denominator = std::pow(6.0, 2.0 / 3.0) * std::pow(z3, 1.0 / 3.0);
    return numerator / denominator;
}
//takes in pressure and density and solves the quartic analytically to get you temperature
Real calcTemperaturePressurePgen(Real rho, Real pres) {
    Real a, mu, mProton, kB, A, B, temp, z1, z2, y;
    mProton = 1.6726 * std::pow(10.0, -24.0);
    kB = 1.3807 * std::pow(10.0, -16.0);
    mu = 4.0 / 3.0;
    a = 7.56 * std::pow(10.0, -15.0);
    A = 3.0 * kB * rho / (a * mu * mProton);
    B = 3.0 * pres / a;
    y = findCubeRoot(A, B);
    temp = std::pow(y, .5) * (std::pow(2.0 * A / std::pow(y * y * y, .5) - 1, .5) - 1) / 2.0;
    return temp;
}

void Mesh::InitUserMeshData(ParameterInput* pin) {
    printf("INITIATING USER MESH DATA\n");
    x1min = pin->GetReal("mesh", "x1min");
    x1max = pin->GetReal("mesh", "x1max");
    nx1 = pin->GetInteger("mesh", "nx1");
    gammagas = pin->GetReal("hydro", "gamma");
    // ambient medium properties
    rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
    p0 = pin->GetOrAddReal("problem", "p0", 1.0);
    vel0 = pin->GetOrAddReal("problem", "vel0", 1.0);
    // gamma law for gas
    gammaGas = pin->GetOrAddReal("hydro", "gamma", 1.66666667);
    // use diode condition for outflow
    diode = pin->GetOrAddBoolean("problem", "diode", false);
    // properties of ejecta
    Rej = pin->GetOrAddReal("problem", "Rej", 1.0);
    Eej = pin->GetOrAddReal("problem", "Eej", 1.0);
    Mej = pin->GetOrAddReal("problem", "Mej", 1.0);
    rhoEj = pin->GetOrAddReal("problem", "rhoEj", 1.0);
    vEj = pin->GetOrAddReal("problem", "vEj", 1.0);
    pressEj = pin->GetOrAddReal("problem", "pressEj", 1.0);
    vmax = pin->GetOrAddReal("problem", "vmax", 1.0);
    opacity = pin->GetOrAddReal("problem", "kappa", .2);
    minRadius = pin->GetOrAddReal("mesh", "x1min", 1.0);
    massFraction = pin->GetOrAddReal("problem", "massFraction", 0.0);
    T0 = pin->GetOrAddReal("radiation", "T_unit", 1.0);
    density0 = pin->GetOrAddReal("radiation", "density_unit", 1.0);
    l0 = pin->GetOrAddReal("radiation", "length_unit", 1.0);
    tempInitial = pin->GetOrAddReal("problem", "T_init", 1.0);
    mu = pin->GetOrAddReal("radiation", "molecular_weight", 4.0/3.0);
    initialTime = pin->GetOrAddReal("problem", "initialTime", 1.0);
    // tell athena to use our user-defined boundaries
    //EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1Wind);
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x1, SNInnerX1Rad);
    //EnrollUserRadBoundaryFunction(BoundaryFace::outer_x1, SNOuterX1Rad);
    EnrollUserExplicitSourceFunction(radioactiveHeating);

    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in octree_test.cpp ProblemGenerator" << std::endl
            << "Spherical-polar coordainates are assumed: " << COORDINATE_SYSTEM << std::endl;
        ATHENA_ERROR(msg);
    }

    return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem generator for testing the octree search
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput* pin) {

    Real r, theta, x, y, dist2, minDist2;
    int index;
    Real kappaES = opacity * density0 * l0;
    Real R = 8.3145 * std::pow(10.0, 7.0);
    Real t0 = Rej / vmax;
    Real v0 = std::pow(T0 * R / mu, .5);
    Real v0sq = 4.0 * Eej / Mej / 3.0;
    Real prefactor = Mej/std::pow(v0sq,1.5)/std::pow(initialTime,3.0)/std::pow(3.14159,1.5);
    int nfreq = pnrrad->nfreq;
    int nang = pnrrad->nang;
    pnrrad->EnrollOpacityFunction(opacityCalc);
        for (int k = ks; k <= ke; k++) {
            for (int j = js; j <= je; j++) {
                for (int i = is; i <= ie; i++) {
                    r = pcoord->x1v(i);
                    theta = pcoord->x2v(j);

                    Real vr   = r/initialTime;
                    Real rho  = prefactor*std::exp(-vr*vr/v0sq);
                    Real pres = 0.7e14*std::pow(rho,1.6666666666667);

                    //x = r * std::cos(theta) * l0;
                    //y = r * std::sin(theta) * l0;

                    Real gasTemp = calcTemperaturePressurePgen(rho, pres);
                    Real currentPres = gasTemp * rho / (1.6036 * std::pow(10.0, -8.0));
                    phydro->u(IDN, k, j, i) = rho;
                    phydro->u(IM1, k, j, i) = rho * vr / v0;
                    phydro->u(IM2, k, j, i) = 0.0;
                    phydro->u(IM3, k, j, i) = 0.0;
                    phydro->u(IEN, k, j, i) = 1.5 * currentPres / ((R / mu) * T0) + 0.5 * rho * (vr / v0) * (vr / v0);
                    Real temperature = 1.6036 * std::pow(10.0, -8.0) * pres / rho;
                    for (int ifr = 0; ifr < pnrrad->nfreq; ++ifr) {
                        pnrrad->sigma_s(k, j, i, ifr)  = rho * kappaES;
                        pnrrad->sigma_pe(k, j, i, ifr) = rho * kappaES;
                        pnrrad->sigma_p(k, j, i, ifr)  = rho * kappaES;
                        for (int n = 0; n < pnrrad->nang; n++) {
                            Real emission = std::pow(gasTemp / T0, 4.0); // (4 * 3.14159265358979);
                            (&(pnrrad->ir(k, j, i, ifr * nang)))[n] = emission;
                        }
                    }
                }
            }
        }
    return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)/
//! \brief Deallocate octree
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput* pin) {
    //printf("Calling octree_final\n");
    //octree_final();
    return;
}

//----------------------------------------------------------------------------------------
//! \fn void SNOuterX1()
//  \brief Sets boundary condition on downstream outflow boundary

void SNOuterX1(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
    Real kappaES, freqScale, temperature, radFactor;
    radFactor = (7.56 * std::pow(10.0, -15.0)) * std::pow(temperature, 4.0) / (3.14159265358979);
    kappaES = opacity * density0 * l0;
    pmb->pnrrad->kappa_es = kappaES;
    bool applyDiode;
    int nfreq = pmb->pnrrad->nfreq;
    Real R = 8.3145 * std::pow(10.0, 7.0);
    Real v0 = std::pow(T0 * R / mu, .5);
    int nang = pmb->pnrrad->nang;
    freqScale = 4.79908742 * std::pow(10.0, -11.0);
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
                temperature = 1.6036 * std::pow(10.0, -8.0) * prim(IEN, k, j, iu) / prim(IDN, k, j, iu);

            }
        }
    }
}

void SNOuterX1Rad(MeshBlock* pmb, Coordinates* pco, NRRadiation* pnrrad, const AthenaArray<Real>& prim, FaceField& b, AthenaArray<Real>& ir,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
    Real kappaES, freqScale, temperature, radFactor;
    radFactor = (7.56 * std::pow(10.0, -15.0)) * std::pow(temperature, 4.0) / (3.14159265358979);
    kappaES = opacity * density0 * l0;
    pnrrad->kappa_es = kappaES;
    bool applyDiode;
    int nfreq = pnrrad->nfreq;
    int nang = pnrrad->nang;
    freqScale = 4.79908742 * std::pow(10.0, -11.0);
    Real R = 8.3145 * std::pow(10.0, 7.0);
    Real v0 = std::pow(T0 * R / mu, .5);
    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {
                temperature = 1.615 * std::pow(10.0, -8.0) * prim(IEN, k, j, iu) * v0 * v0/ prim(IDN, k, j, iu);
                Real emission = std::pow(temperature / T0, 4.0); // (4 * 3.14159265358979);
                for (int ifr = 0; ifr < pnrrad->nfreq; ++ifr) {
                    pnrrad->sigma_s(k, j, i + iu, ifr) = kappaES * prim(IDN, k, j, iu);
                    pnrrad->sigma_a(k, j, i + iu, ifr) = kappaES * prim(IDN, k, j, iu);
                    pnrrad->sigma_pe(k, j, i + iu, ifr) = kappaES * prim(IDN, k, j, iu);
                    pnrrad->sigma_p(k, j, i + iu, ifr) = kappaES * prim(IDN, k, j, iu);
                    for (int n = 0; n < pnrrad->nang; n++) {
                        (&(pnrrad->ir(k, j, i + iu, ifr * nang)))[n] = emission;
                    }
                }

            }
        }
    }
}

Real MyTimeStep(MeshBlock* pmb) {
    Real time = pmb->pmy_mesh->time;
    Real dt = pmb->pmy_mesh->dt;
    Real min_dt = 1e2;
    if (time < .1) {
        min_dt = std::min(dt, 1e-2);
    } // calculate your own time step here
    return min_dt;
}

void SNInnerX1Wind(MeshBlock* pmb, Coordinates* pco, AthenaArray<Real>& prim, FaceField& b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
    int nfreq = pmb->pnrrad->nfreq;
    int nang = pmb->pnrrad->nang;
    Real rho, v, v0sq, pres, prefactor, t0,kappaES, radFactor, freqScale;
    kappaES = opacity * density0 * l0;
    freqScale = 4.79908742 * std::pow(10.0, -11.0);
    pmb->pnrrad->kappa_es = kappaES;
    Real R = 8.3145 * std::pow(10.0, 7.0);
    Real v0 = std::pow(T0 * R / mu, .5);
    if (Rej == 0.0 || Eej == 0.0 || Mej == 0.0 || vmax == 0.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in sn.cpp ProblemGenerator" << std::endl
            << "ejecta parameters not set" << std::endl;
        ATHENA_ERROR(msg);
    }

    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {
                prim(IDN, k, j, il - i) = prim(IDN, k, j, il);
                prim(IM1, k, j, il - i) = prim(IM1, k, j, il);
                prim(IM2, k, j, il - i) = prim(IM2, k, j, il);
                prim(IM3, k, j, il - i) = prim(IM3, k, j, il);
                prim(IEN, k, j, il - i) =  prim(IEN,k,j,il);


            }
        }
    }
}

void SNInnerX1Rad(MeshBlock* pmb, Coordinates* pco, NRRadiation* pnrrad, const AthenaArray<Real>& prim, FaceField& b, AthenaArray<Real>& ir,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
    int nfreq = pnrrad->nfreq;
    int nang = pnrrad->nang;
    Real rho, v, v0sq, pres, prefactor, t0, kappaES, radFactor, freqScale;
    kappaES = opacity * density0 * l0;
    freqScale = 4.79908742 * std::pow(10.0, -11.0);
    pmb->pnrrad->kappa_es = kappaES;
    if (Rej == 0.0 || Eej == 0.0 || Mej == 0.0 || vmax == 0.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in sn.cpp ProblemGenerator" << std::endl
            << "ejecta parameters not set" << std::endl;
        ATHENA_ERROR(msg);
    }

    Real R = 8.3145 * std::pow(10.0, 7.0);
    Real v0 = std::pow(T0 * R / mu, .5);
    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
            for (int i = 1; i <= ngh; ++i) {
                Real temperature = 1.6036 * std::pow(10.0, -8.0) * prim(IEN, k, j, il)*v0*v0 / prim(IDN, k, j, il);
                Real emission = std::pow(temperature / T0, 4.0);//(4*3.14159265358979);
                for (int ifr = 0; ifr < pnrrad->nfreq; ++ifr) {
                    pnrrad->sigma_s(k, j, il - i, ifr) = kappaES * prim(IDN, k, j, il);
                    pnrrad->sigma_a(k, j, il - i, ifr) = kappaES * prim(IDN, k, j, il);
                    pnrrad->sigma_pe(k, j, il - i, ifr) = kappaES * prim(IDN, k, j, il);
                    pnrrad->sigma_p(k, j, il - i, ifr) = kappaES * prim(IDN, k, j, il);
                    for (int n = 0; n < pmb->pnrrad->nang; n++) {
                        (&(pnrrad->ir(k, j, il - i, ifr * nang)))[n] = emission;
                    }
                }

            }
        }
    }
}

void opacityCalc(MeshBlock* pmb, AthenaArray<Real>& prim) {
    int kl = pmb->ks, ku = pmb->ke;
    int jl = pmb->js, ju = pmb->je;
    int il = pmb->is , iu = pmb->ie;
    //if (pmb->block_size.nx2 > 1) {
    //    jl -= 1;
    //    ju += 1;
    //}
    //if (pmb->block_size.nx3 > 1) {
    //    kl -= 1;
    //    ku += 1;
    //}
    Real kappaES = opacity * density0 * l0;
    for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
#pragma omp simd
            for (int i = il; i <= iu; ++i) {
                Real sigma = prim(IDN, k, j, i) * kappaES;
                for (int ifr = 0; ifr < pmb->pnrrad->nfreq; ++ifr) {
                    pmb->pnrrad->sigma_s(k, j, i, ifr) = sigma;
                    pmb->pnrrad->sigma_a(k, j, i, ifr) = sigma;
                    pmb->pnrrad->sigma_pe(k, j, i, ifr) = sigma;
                    pmb->pnrrad->sigma_p(k, j, i, ifr) = sigma;
                }
            }
        }
    }
}

//Implements radioactive heating with Ni-Fe-Co, Ni-Fe has 1.72 MeV
void radioactiveHeating(MeshBlock* pmb, const Real time, const Real dt,
    const AthenaArray<Real>& prim, const AthenaArray<Real>& prim_scalar,
    const AthenaArray<Real>& bcc, AthenaArray<Real>& cons,
    AthenaArray<Real>& cons_scalar) {
    int nfreq = pmb->pnrrad->nfreq;
    int nang = pmb->pnrrad->nang;
    //UPDATE THIS SO IT HAS CORRECT UNITS, JUST EDIT FINAL VALUES AND UPDATE RADIATION LIKE WE HAVE EFFECTIVE TEMPERATURE
    //FOR THIS DO WE JUST ASSUME GAMMA RAYS HEAT GAS AND SO RADIATION FIELD WILL BE UPDATED IN SIM OR DO I HAVE TO UPDATE IT MYSELF?
    //UPDATE GAS ONLY IDC RADIATION
    Real t0 = Rej / vmax;
    Real tau = 757728;
    Real epsilon0 = massFraction * (1.72 * std::pow(10.0, -6.0)) / (tau * 56 * 1.6726 * std::pow(10.0, -24.0)); //epsilon0= X*E_{decay}/(tau*atomicNumber*m_p)
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
            for (int i = pmb->is; i <= pmb->ie; ++i) {
                Real R = 8.3145 * std::pow(10.0, 7.0);
                Real radiationEnergyGrowth = dt*epsilon0 * std::exp(-1.0 * (time + initialTime) / tau) * cons(IDN, k, j, i);
                cons(IEN, k, j, i) += radiationEnergyGrowth / (R * T0/mu);
            }
        }
    }
}
