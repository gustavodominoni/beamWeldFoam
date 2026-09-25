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
    const volScalarField rhoMetal(metalProperty(rho1, rhoB_));
    const volScalarField mu
    (
        limitedAlpha1*rhoMetal*metalNu()
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
        thermalDamper_ =
            2.0*rhoCp/(rhoMetal*metalProperty(cp1_, cpB_) + rho2*cp2_);
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

    // Width of the mushy zone, the inverse of the slope of the liquid
    // fraction with respect to temperature
    const volScalarField deltaTm
    (
        max
        (
            TLiquidus_ - TSolidus_,
            dimensionedScalar(dimTemperature, small)
        )
    );

    const volScalarField rhoByDeltaT(rho/runTime.deltaT());

    // Indicator of the cells changing phase, in which the latent heat
    // source is linearised
    volScalarField phaseChange
    (
        IOobject("phaseChange", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless, 0)
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

        // Evaporative cooling, smoothly switched on over TSmooth around
        // the vapourisation temperature
        const volScalarField evapRamp(vapourRamp(Tvap_));

        const volScalarField QvFull
        (
            evaporativeHeatFlux(Tvap_, Mm_, LatentHeatVap_)
        );

        Qv_ = min(max(evapRamp, scalar(0)), scalar(1))*QvFull;

        // Derivative of the evaporative cooling with respect to
        // temperature, for its linearisation
        tmp<volScalarField> tdQvdT;

        if (metalBVapour_)
        {
            // Sum of the evaporative cooling of the two metals, each
            // weighted by its mole fraction (Raoult's law) and switched on
            // around its own vapourisation temperature
            const volScalarField& xB = metalBMoleFraction_();

            const volScalarField evapRampB(vapourRamp(TvapB_));

            const volScalarField QvFullB
            (
                evaporativeHeatFlux(TvapB_, MmB_, LatentHeatVapB_)
            );

            const volScalarField QvA((1 - xB)*Qv_);
            const volScalarField QvB
            (
                xB*min(max(evapRampB, scalar(0)), scalar(1))*QvFullB
            );

            Qv_ = QvA + QvB;

            if (evaporationLinearisation_)
            {
                tdQvdT =
                    QvA*(LatentHeatVap_*Mm_/(R_*sqr(T_)) - 0.5/T_)
                  + (1 - xB)*pos(evapRamp)*pos(1 - evapRamp)*QvFull/TSmooth_
                  + QvB*(LatentHeatVapB_*MmB_/(R_*sqr(T_)) - 0.5/T_)
                  + xB*pos(evapRampB)*pos(1 - evapRampB)*QvFullB/TSmooth_;
            }
        }
        else if (evaporationLinearisation_)
        {
            tdQvdT =
                Qv_
               *(
                    LatentHeatVap_*Mm_/(R_*sqr(T_))
                  - 0.5/T_
                )
              + pos(evapRamp)*pos(1 - evapRamp)*QvFull/TSmooth_;
        }

        // Temperature consistent with the current liquid fraction
        Tcorr_ = (TLiquidus_ - TSolidus_)*epsilon1_ + TSolidus_;

        fvScalarMatrix TEqn
        (
            TEqnConst
         ==
            fvModels().source(rhoCp, T_)
          - Qv_*magGradAlpha*thermalDamper_
          - TRHS_
        );

        if (latentHeatLinearisation_)
        {
            // Linearise the change of the latent heat source with
            // temperature (Voller and Swaminathan) in the cells which are
            // changing phase, i.e. not fully liquid and above the solidus
            // or not fully solid and below the liquidus. The term vanishes
            // on convergence, where T = Tcorr in all of these cells.
            phaseChange =
                max
                (
                    pos(1 - epsilon1_)*pos(T_ - TSolidus_),
                    pos(epsilon1_)*pos(TLiquidus_ - T_)
                );

            const volScalarField latentCoeff
            (
                phaseChange*LatentHeat_*rhoByDeltaT/deltaTm
            );

            TEqn += fvm::Sp(latentCoeff, T_) - latentCoeff*Tcorr_;
        }

        if (evaporationLinearisation_)
        {
            // Linearise the evaporative cooling, which grows exponentially
            // with temperature, about the current temperature. The term
            // vanishes on convergence and makes the sink implicit, which
            // keeps the temperature bounded for larger time steps.
            const volScalarField dQvdT
            (
                max
                (
                    tdQvdT,
                    dimensionedScalar(Qv_.dimensions()/dimTemperature, 0)
                )
            );

            const volScalarField evapCoeff
            (
                dQvdT*magGradAlpha*thermalDamper_
            );

            TEqn += fvm::Sp(evapCoeff, T_) - evapCoeff*T_;
        }

        TEqn.relax();

        fvConstraints().constrain(TEqn);

        TEqn.solve();

        fvConstraints().constrain(T_);

        // Update the liquid fraction from the temperature, consistently
        // with the linearisation of the latent heat source
        if (latentHeatLinearisation_)
        {
            epsilon1_ =
                max
                (
                    min
                    (
                        epsilon1_
                      + epsilonRel_*phaseChange*(T_ - Tcorr_)/deltaTm,
                        scalar(1)
                    ),
                    scalar(0)
                );
        }
        else
        {
            epsilon1_ =
                max
                (
                    min
                    (
                        epsilon1_
                      + (epsilonRel_*cp_/LatentHeat_)*(T_ - Tcorr_),
                        scalar(1)
                    ),
                    scalar(0)
                );
        }

        residual =
            gMax
            (
                mag
                (
                    epsilon1_.primitiveField()
                  - epsilon1_.prevIter().primitiveField()
                )()
            );

        if (latentHeatLinearisation_)
        {
            // Include the departure of the liquid fraction from its
            // equilibrium value at the new temperature, which catches the
            // cells that start changing phase during this corrector
            residual =
                max
                (
                    residual,
                    gMax
                    (
                        mag
                        (
                            min
                            (
                                max
                                (
                                    (T_ - TSolidus_)/deltaTm,
                                    scalar(0)
                                ),
                                scalar(1)
                            )
                          - epsilon1_
                        )().primitiveField()
                    )
                );
        }
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
