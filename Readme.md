## beamWeldFoam

## Overview
Presented here is the extensible open-source volume-of-fluid (VOF) solver beamWeldFoam, for studying high energy density advanced manufacturing processes. In this implementation the metallic substrate, and shielding gas phase, are treated as in-compressible. The solver fully captures the fusion/melting state transition of the metallic substrate. For the vapourisation of the substrate, the explicit volumetric dilation due to the vapourisation state transition is neglected, instead, a phenomenological recoil pressure term is used to capture the contribution to the momentum and energy fields due to vaporisation events. beamWeldFoam also captures surface tension effects, the temperature dependence of surface tension (Marangoni) effects, latent heat effects due to melting/fusion (and vapourisation), buoyancy effects due to the thermal expansion of the phases using a Boussinesq approximation, momentum damping due to solidification, and a representative heat source description of an incident laser/electron beam heat source. The heat source can also be modified to be representative of arc-welding processes.
The solver is implemented as an OpenFOAM-13 solver module built on the adiabatic two-phase `incompressibleVoF` module (the successor of interFoam) developed by the [OpenFOAM Foundation](https://openfoam.org/). Target applications for beamWeldFoam include:

* Laser Welding
* Electron Beam Welding
* Arc Welding
* Additive Manufacturing

## Installation

The current version of the code is written for the [OpenFOAM-13 libraries](https://openfoam.org/version/13/) from the OpenFOAM Foundation. The code has been developed and tested using an Ubuntu installation, but should work on any operating system capable of installing OpenFOAM. To install the beamWeldFoam solver, first follow the instructions on this page: [OpenFOAM 13 Install](https://openfoam.org/download/13-ubuntu/) to install the OpenFOAM-13 libraries.

Older releases of the code are available for previous OpenFOAM versions: select the `OpenFoam6` branch for OpenFOAM-6 and the `OF10` branch for OpenFOAM-10.

In OpenFOAM-13 beamWeldFoam is implemented as a *solver module* (`solvers::beamWeldFoam`, derived from the standard `incompressibleVoF` module) which is loaded and run by the `foamRun` application. Navigate to a working folder in a shell terminal (with the OpenFOAM-13 environment sourced), clone the git code repository, and build:

```
$ git clone https://github.com/tomflint22/beamWeldFoam.git beamWeldFoam
$ cd beamWeldFoam
$ ./Allwmake
```
This compiles the module into `$FOAM_USER_LIBBIN/libbeamWeldFoam.so` and installs a small `beamWeldFoam` wrapper script into `$FOAM_USER_APPBIN`, which simply calls `foamRun -solver beamWeldFoam`. The installation can be tested using the tutorial cases described below.

## Tutorial cases
Every tutorial contains `Allrun` and `Allclean` scripts. To run a tutorial in serial mode:
```
$ ./Allrun
```
or, step by step:
```
delete any old simulation files, e.g:
$ ./Allclean
Then:
$ cp -r initial 0
$ blockMesh
$ setFields
$ foamRun -solver beamWeldFoam
```
(`beamWeldFoam` can be used in place of `foamRun -solver beamWeldFoam` if the wrapper script installed by `Allwmake` is on your `PATH`; the solver is also selected by the `solver beamWeldFoam;` entry in each `system/controlDict`.)

For parallel deployment, using MPI, following the setFields command:
```
$ decomposePar
$ mpirun -np 6 foamRun -solver beamWeldFoam -parallel >log &
```
for deployment on 6 cores, or simply `./Allrun parallel`. The heat source ray-tracing assumes a structured (blockMesh) mesh and the `simple` or `hierarchical` decomposition methods should be used.

### Case set-up
The solver reads the standard OpenFOAM-13 `incompressibleVoF` case files plus the beamWeldFoam-specific entries:

* `constant/phaseProperties`: the phase names, the surface tension `sigma`, and the interface/vapourisation properties `dsigmadT`, `p0`, `Tvap`, `Mm` and `LatentHeatVap`.
* `constant/physicalProperties.<phase>` (one per phase): the `viscosityModel`, `nu` and `rho` of the phase, and the thermal properties `cp`, `cpsolid`, `kappa`, `kappasolid`, `Tsolidus`, `Tliquidus`, `LatentHeat` and `beta` (`cpsolid` and `kappasolid` default to `cp` and `kappa` if omitted).
* `constant/momentumTransport`: `simulationType laminar;` (or a RAS/LES model).
* `constant/g`: the gravitational acceleration.
* `system/fvSolution`: the `MELTING` sub-dictionary with the liquid-fraction corrector controls and the heat source parameters, in addition to the usual `solvers` and `PIMPLE` entries. The interface compression coefficient is no longer set with `cAlpha` in `fvSolution` but selected with the `div(phi,alpha)` scheme in `system/fvSchemes`, e.g. `Gauss interfaceCompression vanLeer 1;`.
* `initial/` (copied to `0/`): the `alpha.<phase1>`, `p_rgh`, `U` and `Temperature` fields.

The `MELTING` dictionary also accepts the optional switch `writeDiagnostics` (default `true`). When set to `false`, the diagnostic fields (`cp`, `kappa`, `TSolidus`, `TLiquidus`, `LatentHeat`, `beta`, `rhok`, `DC`, `nneps1`, `sourceTerm`, `TRHS`, `ViscousDissipation`, `BeamProfile`, `Num_divU`, `Marangoni`, `pVap`, `Qv`, `ddte1`) are not written. This greatly reduces the output written by large 3D cases. `BeamProfile` is only evaluated at write times.

### Performance tips
The following case settings trade some accuracy for speed. Validate them against a reference run (e.g. the Gallium and Sen & Davies cases) before relying on them:

* `p_rgh`: for large 3D meshes use `GAMG` with `tolerance 1e-8; relTol 0.01;`, keeping `relTol 0` for `p_rghFinal`, instead of `PCG`/`DIC` at `1e-12`.
* `U` and `Temperature`: a tolerance of `1e-8` is normally sufficient; `PBiCGStab` with `DILU` is usually faster than a `symGaussSeidel` smoothSolver for `Temperature`.
* `PIMPLE`: `nCorrectors 3` is usually enough.
* `MELTING`: an `epsilonTolerance` of `1e-6`–`1e-7` usually needs far fewer liquid-fraction corrector iterations than `1e-9`; the number of iterations is reported in the log each time step.
* Output: use `writeFormat binary;`, write less often and set `writeDiagnostics false;`.
* Parallel: use the `simple` or `hierarchical` decomposition with no split in the y (beam) direction.

Cases prepared for the OpenFOAM-6 version of the solver (`constant/transportProperties`, `constant/turbulenceProperties`, `application` entry in `controlDict`, `cAlpha` in `fvSolution`) need to be converted to this layout; the tutorial cases in this repository serve as templates.

### Gallium Melting Case
A commonly used validation case for heat and mass transfer where melting and solidification is involved, is the simulation of Gallium melting in an enclosed container. In this example the beamWeldFoam solver is used to simulate the melting of the Gallium and the subsequent flow due to buoyancy. As time progresses the hot wall on the left-hand-side of the computational domain causes the Gallium in the local vicinity to melt. As the melt volume increases, buoyancy driven flow begins to dominate as the hot liquid Gallium rises and generates vortical flow structures in the liquid. The predicted melt profiles are in excellent agreement with those reported elsewhere, both numerically and experimentally [1].

### Marangoni Flow (Sen and Davies) Case
Another useful validation case for the solver is one in which a 2D cavity is partially filled such that the interface between the phases is initially flat. A temperature gradient is then developed across the domain. This temperature gradient induces a flow tangential to the interface due to the dependence on temperature of the surface tension, aka Marangoni flow. An analytical steady-state solution for the free surface deformation exists for this case [2]. Excellent agreement between the beamWeldFoam solver and the analytical solution is observed.

### Arc Welding Case
In this example a surface heat flux is applied to an Aluminium substrate representative of an arc-welding process. In this scenario, a metallic substrate is present in the domain, between two regions of Argon gas. The heat source is applied at t=0s, and at t=0.25s the power begins to ramp down until at t=0.35s the heat source is fully extinguished. Shortly following the extinction of the heat source the domain fully solidifies. The effect of Marangoni driven flow can clearly be seen in this example, as the surface flows are driven from regions of higher temperature to regions of lower temperature (due to the decrease in surface tension with temperature). Furthermore, once the weld-pool has fully penetrated the domain, surface tension prevents the material from falling out of the bottom of the substrate.

### Beam Welding Case
In this example beamWeldFoam is applied to simulate the power beam welding of a titanium alloy substrate. In this case, Ti6Al4V butt joints welded by a laser beam is simulated and the results are validated with the experimental study [3].

## Algorithm

Initially `foamRun` loads the mesh and constructs the beamWeldFoam solver module, which reads in fields and boundary conditions, reads certain mesh information into arrays (for the heat source application) and selects the momentum transport (turbulence) model (if specified). The main solver loop is then initiated. First, the time step is
dynamically modified to ensure numerical stability. Next, the two-phase fluid mixture properties and turbulence quantities are updated. The discretized phase-fraction equation is then solved for a user-defined number of subtime steps (typically 3) using the multidimensional universal limiter with explicit solution solver [MULES](https://openfoam.org/release/2-3-0/multiphase/). This solver is included in the OpenFOAM library, and performs conservative solution of hyperbolic convective transport equations with defined bounds (0 and 1 for α1). Once the updated phase field is obtained, the program enters the pressure–velocity loop, in which p and u are corrected in an alternating fashion. In this loop T is also solved for, such that he buoyancy predictions are correct for the U and p fields. The process of correcting the pressure and velocity fields in sequence is known as pressure implicit with splitting of operators (PISO). In the OpenFOAM environment, PISO is repeated for multiple iterations at each time step. This process is referred to as merged PISO- semi-implicit method for pressure-linked equations (SIMPLE), or the pressure-velocity loop (PIMPLE) process, where SIMPLE is an iterative pressure–velocity solution algorithm. PIMPLE continues for a user specified number of iterations. 
The main solver loop iterates until program termination. A summary of the simulation algorithm is presented below:
* beamWeldFoam Simulation Algorithm Summary:
  * Initialize simulation data and mesh 
  * WHILE t<t_end DO
  * 1. Update delta_t for stability
  * 2. Phase equation sub-cycle
  * 3. Update interface location for heat source application
  * 4. Update fluid properties
  * 5. PISO Loop
    * 1. Form u equation
    * 2. Energy Transport Loop
      * 1. Solve T equation
      * 2. Update fluid fraction field
      * 3. Re-evaluate source terms due to latent heat
    * 3. PISO
        * 1. Obtain and correct face fluxes
        * 2. Solve p-Poisson equation
        * 3. Correct u
  * 6. Write Fields
  
Two sample tutorial cases, i.e. Gallium Melting, and Sen and Davies cases are in strong agreement with experimental and analytical data available in the literature and serve as the validation cases for the implementation in beamWeldFoam.

## License
OpenFoam, and by extension the beamWeldFoam application, is licensed free and open source only under the [GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). One reason for OpenFOAM’s popularity is that its users are granted the freedom to modify and redistribute the software and have a right of continued free use, within the terms of the GPL.

## Acknowledgements
The work was generously supported by the Engineering and Physical Sciences Research Council (EPSRC) under the ''Cobalt-free Hard-facing for Reactor Systems'' grant EP/T016728/1, and Science Foundation Ireland (SFI), co-funded under European Regional Development Fund and by I-Form industry partners, grant 16/RC/3872.

## Citing This Work
If you use beamWeldFoam in your work. Please use the following to cite our work:

Thomas F. Flint, Gowthaman Parivendhan, Alojz Ivankovic, Michael C. Smith, Philip Cardiff,
beamWeldFoam: Numerical simulation of high energy density fusion and vapourisation-inducing processes,
SoftwareX,
Volume 18,
2022,
101065,
ISSN 2352-7110,
https://doi.org/10.1016/j.softx.2022.101065

## References
* Kay Wittig and Petr A Nikrityuk 2012 IOP Conf. Ser.: Mater. Sci. Eng. 27 012054
* Sen, A., & Davis, S. (1982). Steady thermocapillary flows in two-dimensional slots. Journal of Fluid Mechanics, 121, 163-186. doi:10.1017/S0022112082001840
* Sabina L. Campanelli, Giuseppe Casalino, Michelangelo Mortello, Andrea Angelastro, Antonio Domenico Ludovico, Microstructural Characteristics and Mechanical Properties of Ti6Al4V Alloy Fiber Laser Welds


![visitors](https://visitor-badge.deta.dev/badge?page_id=tomflint22.beamWeldFoam)


