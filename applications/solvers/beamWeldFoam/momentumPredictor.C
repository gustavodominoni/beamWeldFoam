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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::momentumPredictor()
{
    volVectorField& U = U_;

    // Darcy momentum damping in the solid and mushy zone
    DC_ =
        DarcyConstantlarge_*pow((1.0 - epsilon1mask_), 2)
       /(pow(epsilon1mask_, 3) + DarcyConstantsmall_);

    // Marangoni (thermocapillary) force at the interface
    const volVectorField gradAlpha(fvc::grad(alpha1));
    const volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN_));

    Marangoni_ =
        Marangoni_Constant_*(gradT_ - nHatM*(nHatM & gradT_))*mag(gradAlpha);

    if (damperSwitch_)
    {
        // To account for high density difference
        damper_ = 2.0*rho/(mixture.rho1() + mixture.rho2());
    }

    // Recoil pressure
    pVap_ =
        min
        (
            max((T_ - (Tvap_ - (TSmooth_/2.0)))/TSmooth_, scalar(0)),
            scalar(1)
        )
       *0.54*p0_
       *exp(LatentHeatVap_*Mm_*((T_ - Tvap_)/(R_*T_*Tvap_)));

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + divDevTau(U)
      + fvm::Sp(DC_, U)
      - Marangoni_*damper_
     ==
        fvModels().source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    surfaceTensionForce()
                  - buoyancy.ghf*fvc::snGrad(rho*rhok_)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvConstraints().constrain(U);
    }
}


// ************************************************************************* //
