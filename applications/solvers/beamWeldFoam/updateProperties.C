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
#include "fvcAverage.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::updateProperties()
{
    // Mixture specific heat capacity and thermal conductivity blended
    // between the solid and liquid values using the liquid fraction
    cp_ =
        epsilon1_*(alpha1*cp1_ + alpha2*cp2_)
      + (1.0 - epsilon1_)*(alpha1*cp1solid_ + alpha2*cp2solid_);

    kappa_ =
        epsilon1_*(alpha1*kappa1_ + alpha2*kappa2_)
      + (1.0 - epsilon1_)*(alpha1*kappa1solid_ + alpha2*kappa2solid_);

    TSolidus_ = alpha1*Tsolidus1_ + alpha2*Tsolidus2_;
    TLiquidus_ = alpha1*Tliquidus1_ + alpha2*Tliquidus2_;
    LatentHeat_ = alpha1*LatentHeat1_ + alpha2*LatentHeat2_;
    beta_ = alpha1*beta1_ + alpha2*beta2_;

    // Boussinesq density ratio
    rhok_ =
        1.0
      - max
        (
            epsilon1_*beta_*(T_ - TSolidus_),
            dimensionedScalar(dimless, 0)
        );

    // Liquid-fraction interface normal, only needed for output
    if (writeDiagnostics_ && runTime.writeTime())
    {
        const volVectorField gradepsilon1(fvc::grad(epsilon1_));

        nneps1_ = gradepsilon1/(mag(gradepsilon1) + deltaN_);
    }

    // Mask the liquid fraction to remove isolated interface values
    const volScalarField e1temp(fvc::average(epsilon1_));

    scalarField& epsilon1maskI = epsilon1mask_.primitiveFieldRef();
    const scalarField& e1tempI = e1temp.primitiveField();
    const scalarField& epsilon1I = epsilon1_.primitiveField();

    forAll(epsilon1maskI, celli)
    {
        if (e1tempI[celli] <= 0.95)
        {
            epsilon1maskI[celli] = 0.0;
        }
        else
        {
            epsilon1maskI[celli] = epsilon1I[celli];
        }
    }

    epsilon1mask_.correctBoundaryConditions();

    // Buoyancy is only active in the liquid
    buoyancy.gh = epsilon1mask_*(buoyancy.g & mesh.C());
    buoyancy.ghf = fvc::interpolate(epsilon1mask_)*(buoyancy.g & mesh.Cf());
}


// ************************************************************************* //
