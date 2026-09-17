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
#include "meshSearch.H"
#include "DynamicList.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace Foam
{

//- Collect the unique values of the coordinate component cmpt of the
//  cell-centres lying on the line through the first cell-centre along that
//  coordinate direction, i.e. sharing the other two coordinates with the
//  first cell-centre (to within tol)
static scalarList uniqueCoordinates
(
    const vectorField& C,
    const direction cmpt,
    const scalar tol
)
{
    DynamicList<scalar> coords;

    if (C.size())
    {
        const direction cmpt1 = (cmpt + 1) % 3;
        const direction cmpt2 = (cmpt + 2) % 3;

        const vector& C0 = C[0];

        forAll(C, celli)
        {
            const vector& Ci = C[celli];

            if
            (
                mag(Ci[cmpt1] - C0[cmpt1]) < tol
             && mag(Ci[cmpt2] - C0[cmpt2]) < tol
            )
            {
                if (findIndex(coords, Ci[cmpt]) == -1)
                {
                    coords.append(Ci[cmpt]);
                }
            }
        }
    }

    return scalarList(coords);
}

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::findUniqueCoordinates()
{
    const scalar coordTolerance = 1e-7;

    const vectorField& C = mesh.C().primitiveField();

    xlist_ = uniqueCoordinates(C, vector::X, coordTolerance);
    ylist_ = uniqueCoordinates(C, vector::Y, coordTolerance);
    zlist_ = uniqueCoordinates(C, vector::Z, coordTolerance);

    Info<< "Unique cell-centre coordinates: x: " << xlist_.size()
        << " y: " << ylist_.size()
        << " z: " << zlist_.size() << endl;

    if (ylist_.size())
    {
        Info<< "Lowest y co-ordinate: " << min(ylist_) << endl;
    }
}


void Foam::solvers::beamWeldFoam::calcCellGeometry()
{
    const faceList& ff = mesh.faces();
    const pointField& pp = mesh.points();

    forAll(mesh.C(), celli)
    {
        const vector& XYZ = mesh.C()[celli];
        xcoord_[celli] = XYZ.x();
        zcoord_[celli] = XYZ.z();

        const cell& cc = mesh.cells()[celli];
        const labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), Zero);

        forAll(pLabels, pointi)
        {
            pLocal[pointi] = pp[pLabels[pointi]];
        }

        yDim_[celli] =
            Foam::max(pLocal & vector(0, 1, 0))
          - Foam::min(pLocal & vector(0, 1, 0));
    }

    xcoord_.correctBoundaryConditions();
    zcoord_.correctBoundaryConditions();
    yDim_.correctBoundaryConditions();
}


Foam::scalar Foam::solvers::beamWeldFoam::beamIntensity
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar t
) const
{
    using constant::mathematical::pi;

    // Focal position and power, optionally shifted/ramped after tshift
    scalar focnow = focY_;
    scalar Qnow = HS_Q_;

    if (t > tshift_)
    {
        focnow = focY_ - foc_shift_vel_*(t - tshift_);
        Qnow = max(HS_Q_ - Q_ramp_rate_*(t - tshift_), scalar(0));
    }

    // Beam radius diverging away from the focal plane
    const scalar w = HS_a_*Foam::sqrt(1.0 + sqr((y - focnow)/Y_R_));

    // Oscillating beam centre
    const scalar xc =
        HS_bg_
      + Oscillation_Amplitude_*Foam::cos(2.0*pi*Oscillation_Frequency_*t);

    // Travelling beam centre
    const scalar zc = HS_velocity_*t + HS_lg_;

    return
        (2.0*Qnow/(sqr(w)*pi))
       *Foam::exp(-2.0*(sqr((x - xc)/w) + sqr((z - zc)/w)));
}


void Foam::solvers::beamWeldFoam::updateHeatSource()
{
    const scalar t = runTime.value();

    const meshSearch& searchEngine = meshSearch::New(mesh);

    sourceTerm_ = dimensionedScalar(sourceTerm_.dimensions(), 0);

    // Ray-trace along y for each unique (x, z) column of cells to find
    // the first cell containing the substrate and deposit the beam energy
    // in that cell
    forAll(zlist_, zi)
    {
        forAll(xlist_, xi)
        {
            scalar miny = great;

            forAll(ylist_, yi)
            {
                const label celli = searchEngine.findCell
                (
                    point(xlist_[xi], ylist_[yi], zlist_[zi])
                );

                // HS_deposition_cutoff 0.99 for conduction mode,
                // 0.01 for keyhole mode
                if (celli >= 0 && alpha1[celli] > HS_deposition_cutoff_)
                {
                    miny = min(miny, ylist_[yi]);
                }
            }

            if (miny < great)
            {
                const label cellHit = searchEngine.findCell
                (
                    point(xlist_[xi], miny, zlist_[zi])
                );

                if (cellHit >= 0 && yDim_[cellHit] > 1e-12)
                {
                    sourceTerm_[cellHit] =
                        beamIntensity
                        (
                            xcoord_[cellHit],
                            miny,
                            zcoord_[cellHit],
                            t
                        )/yDim_[cellHit];
                }
            }
        }
    }

    sourceTerm_.correctBoundaryConditions();

    // Beam intensity profile for visualisation
    forAll(mesh.C(), celli)
    {
        const vector& XYZ = mesh.C()[celli];
        BeamProfile_[celli] = beamIntensity(XYZ.x(), XYZ.y(), XYZ.z(), t);
    }

    BeamProfile_.correctBoundaryConditions();
}


// ************************************************************************* //
