/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamWeldFoam.H"
#include "fvc.H"
#include "fvm.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::thermophysicalPredictor()
{
    using constant::mathematical::pi;

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    label iter = 0;
    scalar residual = 1;

    TRHS_ = dimensionedScalar(TRHS_.dimensions(), 0);

    const volScalarField rhoCp(rho*cp_);
    const surfaceScalarField rhophicp(fvc::interpolate(cp_)*rhoPhi);

    const volScalarField magGradAlpha(mag(gradAlpha1_));

    // Mixture dynamic viscosity
    const volScalarField limitedAlpha1
    (
        min(max(alpha1, scalar(0)), scalar(1))
    );
    const volScalarField mu
    (
        limitedAlpha1*rho1*mixture.nuModel1().nu()
      + (scalar(1) - limitedAlpha1)*rho2*mixture.nuModel2().nu()
    );

    // Viscous dissipation, which depends only on the velocity and flux and
    // is therefore constant during the liquid-fraction correction
    {
        const volTensorField gradU(fvc::grad(U));
        const volTensorField tau
        (
            mu*gradU + mu*gradU.T() - (2.0/3.0)*mu*fvc::div(phi)*I
        );

        ViscousDissipation_ = tau && gradU;
    }

    if (damperSwitch_)
    {
        thermalDamper_ = 2.0*rhoCp/(rho1*cp1_ + rho2*cp2_);
    }

    // Part of the energy equation which does not change during the
    // liquid-fraction correction
    const fvScalarMatrix TEqnConst
    (
        fvm::ddt(rhoCp, T_)
      + fvm::div(rhophicp, T_)
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T_)
      - fvm::laplacian(kappa_, T_)
      - ViscousDissipation_
     ==
        sourceTerm_
    );

    do
    {
        iter++;

        epsilon1_.storePrevIter();
        T_.storePrevIter();

        // Latent heat source
        TRHS_ =
            LatentHeat_
           *(fvc::ddt(rho, epsilon1_) + fvc::div(rhoPhi, epsilon1_));

        // Evaporative cooling
        Qv_ =
            min
            (
                max((T_ - (Tvap_ - (TSmooth_/2.0)))/TSmooth_, scalar(0)),
                scalar(1)
            )
           *0.82*LatentHeatVap_*Mm_*p0_
           *exp(LatentHeatVap_*Mm_*((T_ - Tvap_)/(R_*T_*Tvap_)))
           /sqrt(2.0*pi*Mm_*R_*T_);

        fvScalarMatrix TEqn
        (
            TEqnConst
         ==
            fvModels().source(rhoCp, T_)
          - Qv_*magGradAlpha*thermalDamper_
          - TRHS_
        );

        TEqn.relax();

        fvConstraints().constrain(TEqn);

        TEqn.solve();

        fvConstraints().constrain(T_);

        // Update the liquid fraction from the temperature
        Tcorr_ = (TLiquidus_ - TSolidus_)*epsilon1_ + TSolidus_;

        epsilon1_ =
            max
            (
                min
                (
                    epsilon1_ + (epsilonRel_*cp_/LatentHeat_)*(T_ - Tcorr_),
                    scalar(1)
                ),
                scalar(0)
            );

        residual =
            gMax
            (
                mag
                (
                    epsilon1_.primitiveField()
                  - epsilon1_.prevIter().primitiveField()
                )()
            );
    }
    while
    (
        (iter < minTCorr_ || residual > epsilonTol_) && iter <= maxTCorr_
    );

    {
        const scalarField depsilon1
        (
            mag
            (
                epsilon1_.primitiveField()
              - epsilon1_.prevIter().primitiveField()
            )
        );

        const scalar meanResidual =
            gSum(depsilon1*mesh.V().primitiveField())
           /gSum(mesh.V().primitiveField());

        Info<< "Correcting epsilon1: " << iter << " iterations"
            << ", mean residual = " << meanResidual
            << ", max residual = " << residual
            << endl;
    }

    ddte1_ = fvc::ddt(epsilon1_);

    T_.correctBoundaryConditions();

    gradT_ = fvc::grad(T_);
}


// ************************************************************************* //
