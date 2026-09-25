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
#include "upwind.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::readMetalB()
{
    typeIOobject<IOdictionary> dictIO
    (
        "physicalProperties.metalB",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!dictIO.headerOk())
    {
        return;
    }

    const IOdictionary dict(dictIO);

    Info<< "Reading the properties of the second metal from "
        << dict.relativeObjectPath() << endl;

    rhoB_.value() = dict.lookup<scalar>("rho");
    nuB_.value() = dict.lookup<scalar>("nu");
    cpB_.value() = dict.lookup<scalar>("cp");
    cpBsolid_.value() = dict.lookupOrDefault<scalar>("cpsolid", cpB_.value());
    kappaB_.value() = dict.lookup<scalar>("kappa");
    kappaBsolid_.value() =
        dict.lookupOrDefault<scalar>("kappasolid", kappaB_.value());
    TsolidusB_.value() = dict.lookup<scalar>("Tsolidus");
    TliquidusB_.value() = dict.lookup<scalar>("Tliquidus");
    LatentHeatB_.value() = dict.lookup<scalar>("LatentHeat");
    betaB_.value() = dict.lookup<scalar>("beta");

    // Vapourisation properties, defaulting to those of metal A
    TvapB_.value() = dict.lookupOrDefault<scalar>("Tvap", Tvap_.value());
    MmB_.value() = dict.lookupOrDefault<scalar>("Mm", Mm_.value());
    LatentHeatVapB_.value() =
        dict.lookupOrDefault<scalar>("LatentHeatVap", LatentHeatVap_.value());

    metalBVapour_ =
        TvapB_.value() != Tvap_.value()
     || MmB_.value() != Mm_.value()
     || LatentHeatVapB_.value() != LatentHeatVap_.value();

    // Surface tension, which can only differ from that of metal A if the
    // latter is constant
    if (dict.found("sigma"))
    {
        if (mixture.isDict("sigma"))
        {
            FatalIOErrorInFunction(dict)
                << "The surface tension of metal B can only be specified "
                << "if that of metal A (sigma in "
                << mixture.relativeObjectPath() << ") is constant"
                << exit(FatalIOError);
        }

        sigmaA_.value() = mixture.lookup<scalar>("sigma");
        sigmaB_.value() = dict.lookup<scalar>("sigma");

        metalBSigma_ = sigmaB_.value() != sigmaA_.value();
    }

    dsigmadTB_.value() =
        dict.lookupOrDefault<scalar>
        (
            "dsigmadT",
            Marangoni_Constant_.value()
        );

    metalBdSigmadT_ = dsigmadTB_.value() != Marangoni_Constant_.value();

    Info<< "    Metal B vapourisation properties "
        << (metalBVapour_ ? "differ from" : "are the same as")
        << " those of metal A" << nl
        << "    Metal B surface tension "
        << (metalBSigma_ ? "differs from" : "is the same as")
        << " that of metal A" << nl
        << "    Metal B surface tension temperature coefficient "
        << (metalBdSigmadT_ ? "differs from" : "is the same as")
        << " that of metal A" << endl;

    alphaMetalB_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "alpha.metalB",
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    metalBFraction_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "metalBFraction",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );

    if (metalBVapour_ || metalBSigma_ || metalBdSigmadT_)
    {
        metalBSurfaceFraction_.reset
        (
            new volScalarField
            (
                IOobject("metalBSurfaceFraction", runTime.name(), mesh),
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
    }

    if (metalBSigma_)
    {
        gradMetalBSurfaceFraction_.reset
        (
            new volVectorField
            (
                IOobject("gradMetalBSurfaceFraction", runTime.name(), mesh),
                mesh,
                dimensionedVector(dimless/dimLength, Zero)
            )
        );
    }

    if (metalBVapour_)
    {
        metalBMoleFraction_.reset
        (
            new volScalarField
            (
                IOobject("metalBMoleFraction", runTime.name(), mesh),
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
    }

    alphaPhiMetalB_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "alphaPhi.metalB",
                runTime.name(),
                mesh
            ),
            alphaPhi1*dimensionedScalar(dimless, 0)
        )
    );

    correctMetalBFraction();
    correctMetalBDensity();
}


void Foam::solvers::beamWeldFoam::correctMetalBFraction()
{
    const volScalarField& alphaB = alphaMetalB_();

    // The fraction is only meaningful where there is metal. Elsewhere it
    // multiplies alpha1, so its value does not matter, but it is kept
    // bounded.
    metalBFraction_() =
        min
        (
            max
            (
                alphaB/max(alpha1, dimensionedScalar(dimless, small)),
                scalar(0)
            ),
            scalar(1)
        );

    if (metalBSurfaceFraction_.valid())
    {
        // The interface cells on the gas side contain little or no metal,
        // so the fraction of metal B there is taken from the average over
        // the neighbouring cells. Without this the surface above metal B
        // would be treated as metal A.
        volScalarField& cS = metalBSurfaceFraction_();

        cS =
            min
            (
                max
                (
                    fvc::average(alphaB)
                   /max
                    (
                        fvc::average(alpha1),
                        dimensionedScalar(dimless, small)
                    ),
                    scalar(0)
                ),
                scalar(1)
            );

        if (gradMetalBSurfaceFraction_.valid())
        {
            gradMetalBSurfaceFraction_() = fvc::grad(cS);
        }

        if (metalBMoleFraction_.valid())
        {
            // Convert the volume fraction to the mole fraction for Raoult's
            // law
            const volScalarField molesB(cS*rhoB_/MmB_);

            metalBMoleFraction_() =
                molesB/((1 - cS)*mixture.rho1()/Mm_ + molesB);
        }
    }
}


void Foam::solvers::beamWeldFoam::transportMetalB()
{
    volScalarField& alphaB = alphaMetalB_();

    // Metal B is carried by the metal flux alphaPhi1 of this time step,
    // in proportion to the upwind fraction of the metal that is metal B,
    // which keeps it consistent with the phase fraction alpha1. Starting
    // from the old-time value allows repeated PIMPLE iterations.
    const volScalarField fractionOld
    (
        min
        (
            max
            (
                alphaB.oldTime()
               /max(alpha1.oldTime(), dimensionedScalar(dimless, small)),
                scalar(0)
            ),
            scalar(1)
        )
    );

    alphaPhiMetalB_() =
        alphaPhi1*upwind<scalar>(mesh, alphaPhi1).interpolate(fractionOld);

    alphaB =
        alphaB.oldTime()
      - runTime.deltaT()*fvc::div(alphaPhiMetalB_());

    // Keep metal B within the metal
    alphaB = min(max(alphaB, scalar(0)), alpha1);
    alphaB.correctBoundaryConditions();

    correctMetalBFraction();
}


void Foam::solvers::beamWeldFoam::correctMetalBDensity()
{
    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    const volScalarField& alphaB = alphaMetalB_();

    // The mixture density and viscosity are owned by the two-phase mixture,
    // which only exposes them as constant, so they are corrected through the
    // object registry
    volScalarField& rhoMix = mesh.lookupObjectRef<volScalarField>("rho");
    volScalarField& nuMix = mesh.lookupObjectRef<volScalarField>("nu");

    rhoMix = alpha1*rho1 + alpha2*rho2 + alphaB*(rhoB_ - rho1);

    // Mass flux consistent with the transport of alpha1 and alpha.metalB
    if (alphaPhiMetalB_.valid())
    {
        rhoPhi =
            alphaPhi1*rho1 + alphaPhi2*rho2
          + alphaPhiMetalB_()*(rhoB_ - rho1);
    }

    // Average kinematic viscosity calculated from the dynamic viscosity, as
    // in incompressibleTwoPhaseVoFMixture
    const volScalarField limitedAlpha1
    (
        min(max(alpha1, scalar(0)), scalar(1))
    );
    const volScalarField rhoMetal(metalProperty(rho1, rhoB_));

    nuMix =
        (
            limitedAlpha1*rhoMetal*metalNu()
          + (scalar(1) - limitedAlpha1)*rho2*mixture.nuModel2().nu()
        )/(limitedAlpha1*rhoMetal + (scalar(1) - limitedAlpha1)*rho2);
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::metalProperty
(
    const dimensionedScalar& a,
    const dimensionedScalar& b
) const
{
    if (metalBFraction_.valid())
    {
        return a + metalBFraction_()*(b - a);
    }
    else
    {
        return volScalarField::New(a.name(), mesh, a);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::surfaceMetalProperty
(
    const dimensionedScalar& a,
    const dimensionedScalar& b
) const
{
    return a + metalBSurfaceFraction_()*(b - a);
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::vapourRamp(const dimensionedScalar& Tvap) const
{
    return (T_ - (Tvap - (TSmooth_/2.0)))/TSmooth_;
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::vapourPressure
(
    const dimensionedScalar& Tvap,
    const dimensionedScalar& Mm,
    const dimensionedScalar& LatentHeatVap
) const
{
    return p0_*exp(LatentHeatVap*Mm*((T_ - Tvap)/(R_*T_*Tvap)));
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::evaporativeHeatFlux
(
    const dimensionedScalar& Tvap,
    const dimensionedScalar& Mm,
    const dimensionedScalar& LatentHeatVap
) const
{
    using constant::mathematical::pi;

    return
        0.82*LatentHeatVap*Mm*p0_
       *exp(LatentHeatVap*Mm*((T_ - Tvap)/(R_*T_*Tvap)))
       /sqrt(2.0*pi*Mm*R_*T_);
}


Foam::tmp<Foam::volScalarField>
Foam::solvers::beamWeldFoam::metalNu() const
{
    if (metalBFraction_.valid())
    {
        const dimensionedScalar& rho1 = mixture.rho1();

        // Blend the dynamic viscosities
        return
            (
                rho1*mixture.nuModel1().nu()
              + metalBFraction_()*(rhoB_*nuB_ - rho1*mixture.nuModel1().nu())
            )/metalProperty(rho1, rhoB_);
    }
    else
    {
        return mixture.nuModel1().nu();
    }
}


// ************************************************************************* //
