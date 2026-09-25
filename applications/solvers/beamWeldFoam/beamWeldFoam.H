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

Class
    Foam::solvers::beamWeldFoam

Description
    Solver module for 2 incompressible, non-isothermal immiscible fluids
    experiencing fusion state changes using a VOF (volume of fluid)
    phase-fraction based interface capturing approach, with optional mesh
    motion and mesh topology changes including adaptive re-meshing.

    Captures the melting and fusion of an alloy, Marangoni (thermocapillary)
    flow, a phenomenological recoil pressure and evaporative cooling model
    for the vapourisation of the substrate, Boussinesq buoyancy, momentum
    damping in the mushy zone and a representative moving (optionally
    oscillating, diverging Gaussian) laser/electron beam heat source.

    The module extends the standard OpenFOAM incompressibleVoF solver module
    and is run with foamRun:

    \verbatim
        foamRun -solver beamWeldFoam
    \endverbatim

    or by setting \c solver beamWeldFoam; in the controlDict.

    Credit:

    Tom Flint, University of Manchester
    Gowthaman Parivendhan, University College Dublin
    Philip Cardiff, University College Dublin

SourceFiles
    beamWeldFoam.C
    updateProperties.C
    heatSource.C
    momentumPredictor.C
    thermophysicalPredictor.C
    pressureCorrector.C

See also
    Foam::solvers::incompressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef beamWeldFoam_H
#define beamWeldFoam_H

#include "incompressibleVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class beamWeldFoam Declaration
\*---------------------------------------------------------------------------*/

class beamWeldFoam
:
    public incompressibleVoF
{

protected:

    // Phase thermal properties
    //  Read from constant/physicalProperties.<phase>

        //- Liquid specific heat capacity of phase 1 and 2
        dimensionedScalar cp1_;
        dimensionedScalar cp2_;

        //- Solid specific heat capacity of phase 1 and 2
        dimensionedScalar cp1solid_;
        dimensionedScalar cp2solid_;

        //- Liquid thermal conductivity of phase 1 and 2
        dimensionedScalar kappa1_;
        dimensionedScalar kappa2_;

        //- Solid thermal conductivity of phase 1 and 2
        dimensionedScalar kappa1solid_;
        dimensionedScalar kappa2solid_;

        //- Solidus temperature of phase 1 and 2
        dimensionedScalar Tsolidus1_;
        dimensionedScalar Tsolidus2_;

        //- Liquidus temperature of phase 1 and 2
        dimensionedScalar Tliquidus1_;
        dimensionedScalar Tliquidus2_;

        //- Latent heat of fusion of phase 1 and 2
        dimensionedScalar LatentHeat1_;
        dimensionedScalar LatentHeat2_;

        //- Thermal expansion coefficient of phase 1 and 2
        dimensionedScalar beta1_;
        dimensionedScalar beta2_;


    // Optional second metal (metal B)
    //  Read from constant/physicalProperties.metalB if present. Phase 1 is
    //  then a mixture of metal A (physicalProperties.<phase1>) and metal B,
    //  the volume fraction of which is alpha.metalB.

        //- Density of metal B
        dimensionedScalar rhoB_;

        //- Kinematic viscosity of liquid metal B
        dimensionedScalar nuB_;

        //- Liquid and solid specific heat capacity of metal B
        dimensionedScalar cpB_;
        dimensionedScalar cpBsolid_;

        //- Liquid and solid thermal conductivity of metal B
        dimensionedScalar kappaB_;
        dimensionedScalar kappaBsolid_;

        //- Solidus and liquidus temperature of metal B
        dimensionedScalar TsolidusB_;
        dimensionedScalar TliquidusB_;

        //- Latent heat of fusion of metal B
        dimensionedScalar LatentHeatB_;

        //- Thermal expansion coefficient of metal B
        dimensionedScalar betaB_;

        //- Vapourisation temperature, molar mass and latent heat of
        //  vapourisation of metal B, defaulting to those of metal A
        dimensionedScalar TvapB_;
        dimensionedScalar MmB_;
        dimensionedScalar LatentHeatVapB_;

        //- Surface tension of metal A and B, only used if they differ
        dimensionedScalar sigmaA_;
        dimensionedScalar sigmaB_;

        //- Temperature coefficient of surface tension of metal B,
        //  defaulting to that of metal A
        dimensionedScalar dsigmadTB_;

        //- Switch for the vapourisation properties of metal B differing
        //  from those of metal A
        bool metalBVapour_;

        //- Switch for the surface tension of metal B differing from that
        //  of metal A
        bool metalBSigma_;

        //- Switch for the temperature coefficient of surface tension of
        //  metal B differing from that of metal A
        bool metalBdSigmadT_;

        //- Volume fraction of metal B
        autoPtr<volScalarField> alphaMetalB_;

        //- Fraction of the metal that is metal B, alpha.metalB/alpha1
        autoPtr<volScalarField> metalBFraction_;

        //- Volumetric flux of metal B
        autoPtr<surfaceScalarField> alphaPhiMetalB_;

        //- Fraction of the metal that is metal B at the interface, which
        //  extends the fraction to the gas side of the interface cells
        autoPtr<volScalarField> metalBSurfaceFraction_;

        //- Gradient of the fraction of metal B at the interface, for the
        //  solutal Marangoni force
        autoPtr<volVectorField> gradMetalBSurfaceFraction_;

        //- Mole fraction of metal B at the interface, for the vapour
        //  pressure (Raoult's law)
        autoPtr<volScalarField> metalBMoleFraction_;


    // Interface and vapourisation properties
    //  Read from constant/phaseProperties

        //- Temperature coefficient of surface tension dsigma/dT
        dimensionedScalar Marangoni_Constant_;

        //- Atmospheric pressure
        dimensionedScalar p0_;

        //- Vapourisation temperature
        dimensionedScalar Tvap_;

        //- Molar mass
        dimensionedScalar Mm_;

        //- Latent heat of vapourisation
        dimensionedScalar LatentHeatVap_;

        //- Universal gas constant
        const dimensionedScalar R_;

        //- Temperature range over which the vapourisation terms are
        //  smoothly switched on
        const dimensionedScalar TSmooth_;

        //- Darcy (mushy-zone) momentum damping constants
        const dimensionedScalar DarcyConstantlarge_;
        const dimensionedScalar DarcyConstantsmall_;

        //- Stabilisation for normalisation of gradients
        const dimensionedScalar deltaN_;


    // Thermophysical fields

        //- Temperature
        volScalarField T_;

        //- Mixture specific heat capacity
        volScalarField cp_;

        //- Mixture thermal conductivity
        volScalarField kappa_;

        //- Mixture solidus temperature
        volScalarField TSolidus_;

        //- Mixture liquidus temperature
        volScalarField TLiquidus_;

        //- Mixture latent heat of fusion
        volScalarField LatentHeat_;

        //- Mixture thermal expansion coefficient
        volScalarField beta_;

        //- Boussinesq density ratio 1 - beta*(T - TSolidus)
        volScalarField rhok_;

        //- Darcy momentum damping coefficient
        volScalarField DC_;

        //- Liquid fraction
        volScalarField epsilon1_;

        //- Liquid fraction masked to remove isolated interface values
        volScalarField epsilon1mask_;

        //- Normalised liquid-fraction gradient
        volVectorField nneps1_;

        //- Temperature gradient
        volVectorField gradT_;

        //- Volumetric heat source from the beam
        volScalarField sourceTerm_;

        //- Latent heat source term
        volScalarField TRHS_;

        //- Temperature consistent with the current liquid fraction
        volScalarField Tcorr_;

        //- Viscous dissipation
        volScalarField ViscousDissipation_;

        //- Beam intensity profile (for visualisation)
        volScalarField BeamProfile_;

        //- Divergence of the volumetric flux (for diagnostics)
        volScalarField Num_divU_;

        //- Marangoni (thermocapillary) force
        volVectorField Marangoni_;

        //- Recoil pressure
        volScalarField pVap_;

        //- Evaporative cooling heat flux
        volScalarField Qv_;

        //- Rate of change of liquid fraction
        volScalarField ddte1_;

        //- Density-ratio damper for the interface forces
        volScalarField damper_;

        //- Density-ratio damper for the interface heat sink
        volScalarField thermalDamper_;

        //- Phase-fraction gradient, cached once per PIMPLE iteration
        volVectorField gradAlpha1_;


    // Cell geometry fields used by the heat source ray-tracing

        //- Cell-centre x-coordinate
        volScalarField xcoord_;

        //- Cell-centre z-coordinate
        volScalarField zcoord_;

        //- Cell extent in the y-direction
        volScalarField yDim_;

        //- Unique cell-centre x, y and z coordinates of the (structured) mesh
        scalarList xlist_;
        scalarList ylist_;
        scalarList zlist_;

        //- Cells of each (x, z) column of the structured mesh, ordered by
        //  increasing cell-centre y-coordinate
        labelListList cellColumns_;


    // Controls
    //  Read from the MELTING sub-dictionary of fvSolution

        //- Minimum and maximum number of temperature/liquid-fraction correctors
        label minTCorr_;
        label maxTCorr_;

        //- Liquid-fraction convergence tolerance
        scalar epsilonTol_;

        //- Liquid-fraction relaxation factor
        scalar epsilonRel_;

        //- Heat source parameters
        scalar HS_a_;
        scalar HS_bg_;
        scalar HS_velocity_;
        scalar HS_lg_;
        scalar HS_Q_;
        scalar HS_deposition_cutoff_;
        scalar Oscillation_Amplitude_;
        scalar Oscillation_Frequency_;
        scalar Y_R_;
        scalar focY_;
        scalar foc_shift_vel_;
        scalar Q_ramp_rate_;
        scalar tshift_;

        //- Switch to damp the interface forces by the density ratio
        bool damperSwitch_;

        //- Switch to write the diagnostic fields
        bool writeDiagnostics_;

        //- Switch to linearise the latent heat source implicitly in the
        //  energy equation
        bool latentHeatLinearisation_;

        //- Switch to linearise the evaporative cooling implicitly in the
        //  energy equation
        bool evaporationLinearisation_;


    // Protected Member Functions

        //- Read controls
        virtual bool read();

        //- Find the unique cell-centre coordinates for the heat source
        //  ray-tracing
        void findUniqueCoordinates(const bool verbose = true);

        //- Calculate the cell geometry fields
        void calcCellGeometry();

        //- Build the ordered cell columns for the heat source ray-tracing
        //  from the unique cell-centre coordinates
        void buildCellColumns();

        //- Set the write option of the diagnostic fields from
        //  writeDiagnostics_
        void setDiagnosticsWriteOpt();

        //- Read the properties and volume fraction of metal B, if present
        void readMetalB();

        //- Update the fraction of the metal that is metal B from the volume
        //  fractions
        void correctMetalBFraction();

        //- Transport the volume fraction of metal B with the metal flux
        void transportMetalB();

        //- Correct the mixture density, viscosity and mass flux for the
        //  density and viscosity of metal B
        void correctMetalBDensity();

        //- Return the property of the metal, blended between the values a
        //  of metal A and b of metal B by the fraction of metal B
        tmp<volScalarField> metalProperty
        (
            const dimensionedScalar& a,
            const dimensionedScalar& b
        ) const;

        //- Return the kinematic viscosity of the liquid metal
        tmp<volScalarField> metalNu() const;

        //- Return the property of the metal at the interface, blended
        //  between the values a of metal A and b of metal B by the fraction
        //  of metal B at the interface
        tmp<volScalarField> surfaceMetalProperty
        (
            const dimensionedScalar& a,
            const dimensionedScalar& b
        ) const;

        //- Return the smooth switch of the vapourisation terms around the
        //  vapourisation temperature Tvap
        tmp<volScalarField> vapourRamp(const dimensionedScalar& Tvap) const;

        //- Return the saturated vapour pressure of a metal
        tmp<volScalarField> vapourPressure
        (
            const dimensionedScalar& Tvap,
            const dimensionedScalar& Mm,
            const dimensionedScalar& LatentHeatVap
        ) const;

        //- Return the evaporative heat flux of a pure metal
        //  (Hertz-Knudsen-Langmuir), before the vapourisation switch
        tmp<volScalarField> evaporativeHeatFlux
        (
            const dimensionedScalar& Tvap,
            const dimensionedScalar& Mm,
            const dimensionedScalar& LatentHeatVap
        ) const;

        //- Update the mixture thermophysical properties
        void updateProperties();

        //- Beam intensity at the given location and time
        scalar beamIntensity
        (
            const scalar x,
            const scalar y,
            const scalar z,
            const scalar t
        ) const;

        //- Update the beam heat source
        void updateHeatSource();

        //- Return the interface surface tension and recoil pressure force
        //  for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;


public:

    //- Runtime type information
    TypeName("beamWeldFoam");


    // Constructors

        //- Construct from region mesh
        beamWeldFoam(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        beamWeldFoam(const beamWeldFoam&) = delete;


    //- Destructor
    virtual ~beamWeldFoam();


    // Member Functions

        //- Corrections that follow mesh motion
        //  Rebuilds the cached ray-tracing data, which mesh motion
        //  invalidates
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation and liquid fraction
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const beamWeldFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
