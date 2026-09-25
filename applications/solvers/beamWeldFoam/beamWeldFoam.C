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
#include "findRefCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(beamWeldFoam, 0);
    addToRunTimeSelectionTable(solver, beamWeldFoam, fvMesh);
}
}


// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace Foam
{

//- Construct an IOobject for a field of the given name in the current time
static IOobject fieldIO
(
    const word& name,
    const fvMesh& mesh,
    const IOobject::readOption r,
    const IOobject::writeOption w
)
{
    return IOobject(name, mesh.time().name(), mesh, r, w);
}


//- Read a dimensioned property from a dictionary
static dimensionedScalar readProperty
(
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict,
    const word& key
)
{
    return dimensionedScalar(name, dims, dict.lookup<scalar>(key));
}


//- Read a dimensioned property from a dictionary, falling back to the value
//  of another key if the requested key is not present
static dimensionedScalar readPropertyOrDefault
(
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict,
    const word& key,
    const word& defaultKey
)
{
    return dimensionedScalar
    (
        name,
        dims,
        dict.lookupOrDefault<scalar>(key, dict.lookup<scalar>(defaultKey))
    );
}

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::solvers::beamWeldFoam::read()
{
    incompressibleVoF::read();

    const dictionary& meltingDict = mesh.solution().subDict("MELTING");

    minTCorr_ = meltingDict.lookup<label>("minTempCorrector");
    maxTCorr_ = meltingDict.lookup<label>("maxTempCorrector");

    epsilonTol_ = meltingDict.lookup<scalar>("epsilonTolerance");
    epsilonRel_ = meltingDict.lookup<scalar>("epsilonRelaxation");

    HS_a_ = meltingDict.lookup<scalar>("HS_a");
    HS_bg_ = meltingDict.lookup<scalar>("HS_bg");
    HS_velocity_ = meltingDict.lookup<scalar>("HS_velocity");
    HS_lg_ = meltingDict.lookup<scalar>("HS_lg");
    HS_Q_ = meltingDict.lookup<scalar>("HS_Q");

    HS_deposition_cutoff_ = meltingDict.lookup<scalar>("HS_deposition_cutoff");
    Oscillation_Amplitude_ = meltingDict.lookup<scalar>("Oscillation_Amplitude");
    Oscillation_Frequency_ = meltingDict.lookup<scalar>("Oscillation_Frequency");

    Y_R_ = meltingDict.lookup<scalar>("Y_R");
    focY_ = meltingDict.lookup<scalar>("focY");
    foc_shift_vel_ = meltingDict.lookup<scalar>("foc_shift_vel");
    Q_ramp_rate_ = meltingDict.lookup<scalar>("Q_ramp_rate");
    tshift_ = meltingDict.lookup<scalar>("tshift");

    damperSwitch_ = meltingDict.lookupOrDefault<bool>("damperSwitch", false);

    writeDiagnostics_ =
        meltingDict.lookupOrDefault<bool>("writeDiagnostics", true);

    latentHeatLinearisation_ =
        meltingDict.lookupOrDefault<bool>("latentHeatLinearisation", true);

    evaporationLinearisation_ =
        meltingDict.lookupOrDefault<bool>("evaporationLinearisation", true);

    setDiagnosticsWriteOpt();

    return true;
}


void Foam::solvers::beamWeldFoam::setDiagnosticsWriteOpt()
{
    const IOobject::writeOption w =
        writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;

    // gradT is always written: the momentum predictor uses the value from
    // the previous time step, so it is needed to restart exactly
    cp_.writeOpt() = w;
    kappa_.writeOpt() = w;
    TSolidus_.writeOpt() = w;
    TLiquidus_.writeOpt() = w;
    LatentHeat_.writeOpt() = w;
    beta_.writeOpt() = w;
    rhok_.writeOpt() = w;
    DC_.writeOpt() = w;
    nneps1_.writeOpt() = w;
    sourceTerm_.writeOpt() = w;
    TRHS_.writeOpt() = w;
    ViscousDissipation_.writeOpt() = w;
    BeamProfile_.writeOpt() = w;
    Num_divU_.writeOpt() = w;
    Marangoni_.writeOpt() = w;
    pVap_.writeOpt() = w;
    Qv_.writeOpt() = w;
    ddte1_.writeOpt() = w;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::beamWeldFoam::beamWeldFoam(fvMesh& mesh)
:
    incompressibleVoF(mesh),

    cp1_
    (
        readProperty("cp1", dimSpecificHeatCapacity, mixture.nuModel1(), "cp")
    ),
    cp2_
    (
        readProperty("cp2", dimSpecificHeatCapacity, mixture.nuModel2(), "cp")
    ),

    cp1solid_
    (
        readPropertyOrDefault
        (
            "cp1solid",
            dimSpecificHeatCapacity,
            mixture.nuModel1(),
            "cpsolid",
            "cp"
        )
    ),
    cp2solid_
    (
        readPropertyOrDefault
        (
            "cp2solid",
            dimSpecificHeatCapacity,
            mixture.nuModel2(),
            "cpsolid",
            "cp"
        )
    ),

    kappa1_
    (
        readProperty
        (
            "kappa1",
            dimThermalConductivity,
            mixture.nuModel1(),
            "kappa"
        )
    ),
    kappa2_
    (
        readProperty
        (
            "kappa2",
            dimThermalConductivity,
            mixture.nuModel2(),
            "kappa"
        )
    ),

    kappa1solid_
    (
        readPropertyOrDefault
        (
            "kappa1solid",
            dimThermalConductivity,
            mixture.nuModel1(),
            "kappasolid",
            "kappa"
        )
    ),
    kappa2solid_
    (
        readPropertyOrDefault
        (
            "kappa2solid",
            dimThermalConductivity,
            mixture.nuModel2(),
            "kappasolid",
            "kappa"
        )
    ),

    Tsolidus1_
    (
        readProperty("Tsolidus1", dimTemperature, mixture.nuModel1(), "Tsolidus")
    ),
    Tsolidus2_
    (
        readProperty("Tsolidus2", dimTemperature, mixture.nuModel2(), "Tsolidus")
    ),

    Tliquidus1_
    (
        readProperty
        (
            "Tliquidus1",
            dimTemperature,
            mixture.nuModel1(),
            "Tliquidus"
        )
    ),
    Tliquidus2_
    (
        readProperty
        (
            "Tliquidus2",
            dimTemperature,
            mixture.nuModel2(),
            "Tliquidus"
        )
    ),

    LatentHeat1_
    (
        readProperty
        (
            "LatentHeat1",
            dimEnergy/dimMass,
            mixture.nuModel1(),
            "LatentHeat"
        )
    ),
    LatentHeat2_
    (
        readProperty
        (
            "LatentHeat2",
            dimEnergy/dimMass,
            mixture.nuModel2(),
            "LatentHeat"
        )
    ),

    beta1_
    (
        readProperty("beta1", dimless/dimTemperature, mixture.nuModel1(), "beta")
    ),
    beta2_
    (
        readProperty("beta2", dimless/dimTemperature, mixture.nuModel2(), "beta")
    ),

    rhoB_("rhoB", dimDensity, 0),
    nuB_("nuB", dimKinematicViscosity, 0),
    cpB_("cpB", dimSpecificHeatCapacity, 0),
    cpBsolid_("cpBsolid", dimSpecificHeatCapacity, 0),
    kappaB_("kappaB", dimThermalConductivity, 0),
    kappaBsolid_("kappaBsolid", dimThermalConductivity, 0),
    TsolidusB_("TsolidusB", dimTemperature, 0),
    TliquidusB_("TliquidusB", dimTemperature, 0),
    LatentHeatB_("LatentHeatB", dimEnergy/dimMass, 0),
    betaB_("betaB", dimless/dimTemperature, 0),

    Marangoni_Constant_
    (
        readProperty
        (
            "Marangoni_Constant",
            dimensionSet(1, 0, -2, -1, 0),
            mixture,
            "dsigmadT"
        )
    ),
    p0_(readProperty("p0", dimPressure, mixture, "p0")),
    Tvap_(readProperty("Tvap", dimTemperature, mixture, "Tvap")),
    Mm_(readProperty("Mm", dimMass/dimMoles, mixture, "Mm")),
    LatentHeatVap_
    (
        readProperty("LatentHeatVap", dimEnergy/dimMass, mixture, "LatentHeatVap")
    ),

    R_("R", dimensionSet(1, 2, -2, -1, -1), 8.314),
    TSmooth_("TSmooth", dimTemperature, 50.0),
    DarcyConstantlarge_
    (
        "DarcyConstantlarge",
        dimDensity/dimTime,
        scalar(1.0e6)
    ),
    DarcyConstantsmall_("DarcyConstantsmall", dimless, scalar(1.0e-12)),
    deltaN_("deltaN", 1e-8/pow(average(mesh.V()), 1.0/3.0)),

    T_
    (
        fieldIO("Temperature", mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    ),
    cp_
    (
        fieldIO("cp", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimSpecificHeatCapacity, 0)
    ),
    kappa_
    (
        fieldIO("kappa", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimThermalConductivity, 0)
    ),
    TSolidus_
    (
        fieldIO("TSolidus", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimTemperature, 0)
    ),
    TLiquidus_
    (
        fieldIO("TLiquidus", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimTemperature, 1)
    ),
    LatentHeat_
    (
        fieldIO("LatentHeat", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0)
    ),
    beta_
    (
        fieldIO("beta", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTemperature, 0)
    ),
    rhok_
    (
        fieldIO("rhok", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    DC_
    (
        fieldIO("DC", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 1.0e14)
    ),
    epsilon1_
    (
        fieldIO("epsilon1", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    epsilon1mask_
    (
        fieldIO("epsilon1mask", mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    nneps1_
    (
        fieldIO("nneps1", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector(dimless, Zero)
    ),
    gradT_
    (
        fieldIO("gradT", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector(dimTemperature/dimLength, Zero)
    ),
    sourceTerm_
    (
        fieldIO("sourceTerm", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),
    TRHS_
    (
        fieldIO("TRHS", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),
    Tcorr_
    (
        fieldIO("Tcorr", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimTemperature, 0)
    ),
    ViscousDissipation_
    (
        fieldIO
        (
            "ViscousDissipation",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),
    BeamProfile_
    (
        fieldIO("BeamProfile", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),
    Num_divU_
    (
        fieldIO("Num_divU", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTime, 0)
    ),
    Marangoni_
    (
        fieldIO("Marangoni", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector(dimPressure/dimLength, Zero)
    ),
    pVap_
    (
        fieldIO("pVap", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimPressure, 0)
    ),
    Qv_
    (
        fieldIO("Qv", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimPower/dimArea, 0)
    ),
    ddte1_
    (
        fieldIO("ddte1", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    damper_
    (
        fieldIO("damper", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    thermalDamper_
    (
        fieldIO("thermalDamper", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimless, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    gradAlpha1_
    (
        fieldIO("gradAlpha1", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        fvc::grad(alpha1)
    ),

    xcoord_
    (
        fieldIO("xcoord", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimLength, 0)
    ),
    zcoord_
    (
        fieldIO("zcoord", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimLength, 0)
    ),
    yDim_
    (
        fieldIO("yDim", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimLength, 1)
    ),

    minTCorr_(1),
    maxTCorr_(1),
    epsilonTol_(0),
    epsilonRel_(1),
    HS_a_(0),
    HS_bg_(0),
    HS_velocity_(0),
    HS_lg_(0),
    HS_Q_(0),
    HS_deposition_cutoff_(0),
    Oscillation_Amplitude_(0),
    Oscillation_Frequency_(0),
    Y_R_(1),
    focY_(0),
    foc_shift_vel_(0),
    Q_ramp_rate_(0),
    tshift_(0),
    damperSwitch_(false),
    writeDiagnostics_(true),
    latentHeatLinearisation_(true),
    evaporationLinearisation_(true)
{
    // Read the MELTING controls
    read();

    Info<< "\nFinding the unique cell-centre coordinates "
           "for the heat source ray-tracing\n" << endl;
    findUniqueCoordinates();
    buildCellColumns();
    calcCellGeometry();

    // Read the optional second metal and correct the mixture density and
    // viscosity for it
    readMetalB();

    // Initialise the mixture properties and liquid fraction from the
    // initial temperature field
    TSolidus_ =
        alpha1*metalProperty(Tsolidus1_, TsolidusB_) + alpha2*Tsolidus2_;
    TLiquidus_ =
        alpha1*metalProperty(Tliquidus1_, TliquidusB_) + alpha2*Tliquidus2_;

    epsilon1_ =
        max
        (
            min((T_ - TSolidus_)/(TLiquidus_ - TSolidus_), scalar(1)),
            scalar(0)
        );

    updateProperties();

    // Re-evaluate the pressure with the liquid-fraction masked gh
    p = p_rgh + rho*buoyancy.gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureReference().refValue()
          - getRefCellValue(p, pressureReference().refCell())
        );
        p_rgh_ = p - rho*buoyancy.gh;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::beamWeldFoam::~beamWeldFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::beamWeldFoam::motionCorrector()
{
    incompressibleVoF::motionCorrector();

    // The heat-source ray tracing caches the cell-centre coordinates, the
    // cell heights and the lists of unique cell-centre coordinates. All of
    // these describe the mesh geometry, so mesh motion invalidates them and
    // the beam would otherwise continue to be deposited at the old location.
    if (mesh.changing())
    {
        // The ray tracing walks columns of cells along the y-axis, which
        // only identifies a surface cell on a structured mesh of fixed
        // topology. Refinement, layer addition and stitching break that
        // assumption rather than merely moving the data.
        if (mesh.topoChanged())
        {
            FatalErrorInFunction
                << "The beamWeldFoam heat-source ray tracing requires a "
                   "structured mesh of fixed topology, but the mesh topology "
                   "changed." << nl
                << "Mesh topology changes (refinement, layer addition, "
                   "stitching) are not supported by the heat source."
                << exit(FatalError);
        }

        findUniqueCoordinates(false);
        buildCellColumns();
        calcCellGeometry();
    }
}


void Foam::solvers::beamWeldFoam::prePredictor()
{
    // Solve the phase-fraction equation and update the mixture properties
    // and mass flux
    incompressibleVoF::prePredictor();

    // Transport the second metal and correct the mixture density,
    // viscosity and mass flux for it
    if (alphaMetalB_.valid())
    {
        transportMetalB();
        correctMetalBDensity();
    }

    // Cache the phase-fraction gradient used by the momentum and energy
    // equations of this PIMPLE iteration
    gradAlpha1_ = fvc::grad(alpha1);

    // Update the thermophysical properties of the mixture
    updateProperties();

    // Update the beam heat source
    updateHeatSource();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::beamWeldFoam::surfaceTensionForce() const
{
    return
        (
            interface.surfaceTensionForce()
          + fvc::interpolate(pVap_)*fvc::snGrad(alpha1)
        )*fvc::interpolate(damper_);
}


// ************************************************************************* //
