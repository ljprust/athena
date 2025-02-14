//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file boos_input.cpp
//! \brief Problem generator for reading in data from Boos+24

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
#include <unistd.h>   
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
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
//#include "../units/units.hpp"

std::vector<Real> vx_in, vz_in, rho_in, temp_in;
std::vector<Real> ejecta_in, he_in, n_in, o_in, si_in, fe_in;

namespace {
    int NumToRead, nx1;
    Real x1max, deltax;
    Real a, mu, kB, mProton;
    Real rho0, temp0, gammaGas, initialTime;
} // namespace
/*
//gets pressure given temperature and density
Real pressureEqPgen(Real temp, Real rho) {
    return ((1.0 / 3.0) * a * std::pow(temp, 4.0)) + rho * kB * temp / (mu * mProton);
}

//gets internal energy given temperature and density
Real energyEqPgen(Real temp, Real rho) {
    return (a * std::pow(temp, 4.0)) + 1.5 * rho * kB * temp / (mu * mProton);
}

Real gasEnergyEq(Real rho, Real temp) {
    return 1.5 * rho * kB * temp / (mu * mProton);
}

//Helper functions to calculate temperature:
Real findCubeRoot(Real A, Real B) {
    Real z3 = std::pow(81.0 * std::pow(A, 4.0) + 768.0 * std::pow(B, 3.0), 0.5) + 9.0 * A * A;
    Real numerator = std::pow(2.0, 1.0 / 3.0) * std::pow(z3, 2.0 / 3.0) - 8.0 * std::pow(3.0, 1.0 / 3.0) * B;
    Real denominator = std::pow(6.0, 2.0 / 3.0) * std::pow(z3, 1.0 / 3.0);
    return numerator / denominator;
}

//takes in pressure and density and solves the quartic analytically to get you temperature
Real calcTemperaturePressurePgen(Real rho, Real pres) {
    Real A, B, temp, y;
    A = 3.0 * kB * rho / (a * mu * mProton);
    B = 3.0 * pres / a;
    y = findCubeRoot(A, B);
    temp = std::pow(y, 0.5) * (std::pow(2.0 * A / std::pow(y * y * y, 0.5) - 1.0, 0.5) - 1.0) / 2.0;
    return temp;
}
*/
void Mesh::InitUserMeshData(ParameterInput* pin) {
    printf("Starting InitUserMeshData\n");

    // constants
    a  = 7.5646e-15;
    mu = 2.0;
    //R  = 8.314e7;
    mProton = 1.6726e-24;
    kB      = 1.3807e-16;

    NumToRead   = pin->GetInteger("problem", "NumToRead");
    x1max       = pin->GetReal("mesh", "x1max");
    initialTime = pin->GetOrAddReal("problem", "initialTime", 1.0);
    nx1         = pin->GetInteger("mesh", "nx1");

    // ambient medium properties
    rho0  = pin->GetOrAddReal("problem", "rho0",  1.0);
    temp0 = pin->GetOrAddReal("problem", "temp0", 1.0);

    deltax = x1max / (static_cast<double>(nx1));

    if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") != 0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in boos_input.cpp ProblemGenerator" << std::endl
            << "Cylindrical-polar coordainates are assumed: " << COORDINATE_SYSTEM << std::endl;
        ATHENA_ERROR(msg);
    }

    char vxFile[256];
    char vzFile[256];
    char rhoFile[256];
    char tempFile[256];
    char ejectaFile[256];
    char heFile[256];
    char nFile[256];
    char oFile[256];
    char siFile[256];
    char feFile[256];

    sprintf(vxFile,     "athenainput_vx.txt");
    sprintf(vzFile,     "athenainput_vz.txt");
    sprintf(rhoFile,    "athenainput_rho.txt");
    sprintf(tempFile,   "athenainput_temp.txt");
    sprintf(ejectaFile, "athenainput_ejecta.txt");
    sprintf(heFile,     "athenainput_he.txt");
    sprintf(nFile,      "athenainput_n.txt");
    sprintf(oFile,      "athenainput_o.txt");
    sprintf(siFile,     "athenainput_si.txt");
    sprintf(feFile,     "athenainput_fe.txt");

    printf("Opening data files with state variables...\n");
    std::ifstream vxFileRead, vzFileRead, rhoFileRead, tempFileRead;
    std::ifstream ejectaFileRead, heFileRead, nFileRead, oFileRead, siFileRead, feFileRead;
    vxFileRead.open(vxFile);
    vzFileRead.open(vzFile);
    rhoFileRead.open(rhoFile);
    tempFileRead.open(tempFile);
    ejectaFileRead.open(ejectaFile);
    heFileRead.open(heFile);
    nFileRead.open(nFile);
    oFileRead.open(oFile);
    siFileRead.open(siFile);
    feFileRead.open(feFile);

    Real vx, vz, rho, temp;
    Real ejecta_frac, he_frac, n_frac, o_frac, si_frac, fe_frac;
    for (int l = 0; l < NumToRead; l++) {
        vxFileRead     >> vx;
        vzFileRead     >> vz;
        rhoFileRead    >> rho;
        tempFileRead   >> temp;
        ejectaFileRead >> ejecta_frac;
        heFileRead     >> he_frac;
        nFileRead      >> n_frac;
        oFileRead      >> o_frac;
        siFileRead     >> si_frac;
        feFileRead     >> fe_frac;

        vx_in.push_back(vx);
        vz_in.push_back(vz);
        rho_in.push_back(rho);
        temp_in.push_back(temp);
        ejecta_in.push_back(ejecta_frac);
        he_in.push_back(he_frac);
        n_in.push_back(n_frac);
        o_in.push_back(o_frac);
        si_in.push_back(si_frac);
        fe_in.push_back(fe_frac);
    }
    printf("Done reading, closing data files\n");
    vxFileRead.close();
    vzFileRead.close();
    rhoFileRead.close();
    tempFileRead.close();
    ejectaFileRead.close();
    heFileRead.close();
    nFileRead.close();
    oFileRead.close();
    siFileRead.close();
    feFileRead.close();

    return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem generator for testing the octree search
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput* pin) {

    Real dist2, minDist2, vx, vz;
    Real Egas, Erad, Ekin;
    Real Pgas, Prad, gammaGas, beta;
    int index;
    bool isEjecta;

    for (int k = ks; k <= ke; k++) {
        for (int j = js; j <= je; j++) {
            for (int i = is; i <= ie; i++) {
                vx = pcoord->x1v(i)/initialTime;
                vz = pcoord->x3v(k)/initialTime;

                index = -1;
                minDist2 = 4.0 * x1max * x1max;
                for (int l = 0; l < NumToRead; l++) {
                    dist2 = (vx-vx_in[l])*(vx-vx_in[l]) + (vz-vz_in[l])*(vz-vz_in[l]);
                    if (dist2 < minDist2) {
                        minDist2 = dist2;
                        index = l;
                    }
                }
                //if(vx<5.0e9 && vz<5.0e9 && vz>-5.0e9) printf("for vx vz %5.3e %5.3e found neighbor id %d with vx vz %5.3e %5.3e\n",vx,vz,index,vx_in[index],vz_in[index]);

                isEjecta = ejecta_in[index] > 0.5 && minDist2 < 2.0*deltax;

                if ( isEjecta ) {
                    Pgas = rho_in[index]*temp_in[index]*kB/mProton/mu;
                    Prad = a*temp_in[index]*temp_in[index]*temp_in[index]*temp_in[index]/3.0;
                    beta = Pgas/(Pgas+Prad);
                    gammaGas = (32.0-24.0*beta-3.0*beta*beta)/(24.0-21.0*beta);

                    Egas = 1.0/(1.0-gammaGas)*Pgas;
                    Erad = Prad*3.0;
                    Ekin = 0.5*rho_in[index]*(vx_in[index]*vx_in[index]+vz_in[index]*vz_in[index]);

                    phydro->u(IDN,k,j,i) = rho_in[index];
                    phydro->u(IM1,k,j,i) = rho_in[index] * vx_in[index];
                    phydro->u(IM2,k,j,i) = 0.0;
                    phydro->u(IM3,k,j,i) = rho_in[index] * vz_in[index];
                    phydro->u(IEN,k,j,i) = Egas + Erad + Ekin;
                    pscalars->s(0,k,j,i) = ejecta_in[index]; // ejecta fraction
                    pscalars->s(1,k,j,i) = he_in[index] * rho_in[index];
                    pscalars->s(2,k,j,i) = n_in[index]  * rho_in[index];
                    pscalars->s(3,k,j,i) = o_in[index]  * rho_in[index];
                    pscalars->s(4,k,j,i) = si_in[index] * rho_in[index];
                    pscalars->s(5,k,j,i) = fe_in[index] * rho_in[index];
                } else { // ambient medium
                    Egas = 1.5*rho0*temp0*kB/mProton/mu;
                    Erad = a*temp0*temp0*temp0*temp0;

                    phydro->u(IDN,k,j,i) = rho0;
                    phydro->u(IM1,k,j,i) = 0.0;
                    phydro->u(IM2,k,j,i) = 0.0;
                    phydro->u(IM3,k,j,i) = 0.0;
                    phydro->u(IEN,k,j,i) = Egas + Erad;
                    pscalars->s(0,k,j,i) = 0.0;
                    pscalars->s(1,k,j,i) = 0.0;
                    pscalars->s(2,k,j,i) = 0.0;
                    pscalars->s(3,k,j,i) = 0.0;
                    pscalars->s(4,k,j,i) = 0.0;
                    pscalars->s(5,k,j,i) = 0.0;
                }
            }
        }
    }
    return;
}

