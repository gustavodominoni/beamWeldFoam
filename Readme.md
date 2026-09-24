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
for deployment on 6 cores, or simply `./Allrun parallel`. The heat source ray-tracing assumes a structured (blockMesh) mesh. Beam columns split between processors are handled by a parallel reduction, so any decomposition of such a mesh can be used.

### Case set-up
The solver reads the standard OpenFOAM-13 `incompressibleVoF` case files plus the beamWeldFoam-specific entries:

* `constant/phaseProperties`: the phase names, the surface tension `sigma`, and the interface/vapourisation properties `dsigmadT`, `p0`, `Tvap`, `Mm` and `LatentHeatVap`.
* `constant/physicalProperties.<phase>` (one per phase): the `viscosityModel`, `nu` and `rho` of the phase, and the thermal properties `cp`, `cpsolid`, `kappa`, `kappasolid`, `Tsolidus`, `Tliquidus`, `LatentHeat` and `beta` (`cpsolid` and `kappasolid` default to `cp` and `kappa` if omitted).
* `constant/momentumTransport`: `simulationType laminar;` (or a RAS/LES model).
* `constant/g`: the gravitational acceleration.
* `system/fvSolution`: the `MELTING` sub-dictionary with the liquid-fraction corrector controls and the heat source parameters, in addition to the usual `solvers` and `PIMPLE` entries. The interface compression coefficient is no longer set with `cAlpha` in `fvSolution` but selected with the `div(phi,alpha)` scheme in `system/fvSchemes`, e.g. `Gauss interfaceCompression vanLeer 1;`.
* `initial/` (copied to `0/`): the `alpha.<phase1>`, `p_rgh`, `U` and `Temperature` fields.

The `MELTING` dictionary also accepts the optional switch `writeDiagnostics` (default `true`). When set to `false`, the diagnostic fields (`cp`, `kappa`, `TSolidus`, `TLiquidus`, `LatentHeat`, `beta`, `rhok`, `DC`, `nneps1`, `sourceTerm`, `TRHS`, `ViscousDissipation`, `BeamProfile`, `Num_divU`, `Marangoni`, `pVap`, `Qv`, `ddte1`) are not written. This greatly reduces the output written by large 3D cases. `BeamProfile` is only evaluated at write times.

The liquid-fraction corrector linearises the latent heat source implicitly in the energy equation in the cells that are changing phase (the Voller–Swaminathan source-based method), and updates the liquid fraction consistently with this linearisation. This converges to the same solution as the explicit update of earlier versions, but in far fewer corrector iterations; convergence is judged on both the change in the liquid fraction and its departure from the equilibrium value at the new temperature. The earlier explicit scheme can be restored with `latentHeatLinearisation false;` in `MELTING`.

The evaporative cooling, which grows exponentially with temperature, is likewise linearised implicitly about the current temperature in the energy equation (`evaporationLinearisation`, default `true`). This leaves the solution unchanged until evaporation starts, but keeps the temperature bounded once it does: with the explicit treatment the EB_3D case failed with a negative temperature shortly after the onset of vapourisation (t = 6.3e-5 s), even at its original time step of 2.5e-8 s.

### Performance tips
The following case settings trade some accuracy for speed. Validate them against a reference run (e.g. the Gallium and Sen & Davies cases) before relying on them:

* `p_rgh`: a plain `GAMG` solver with the `GaussSeidel` smoother can stall on the first time step and, as in the ArcCase tutorial at t ≈ 0.0113 s, diverge. `PCG` preconditioned by `GAMG` with the `DICGaussSeidel` smoother (as now used by the ArcCase and PowderBed2D tutorials) is robust and needs at most a few tens of iterations. Switching PowderBed3D from `PCG`/`DIC` to `GAMG` gave no gain.
* `Temperature`: `PBiCGStab` with `DILU` converges in far fewer iterations than a `symGaussSeidel` smoothSolver at the same tolerance (PowderBed3D: 1.3× faster per time step with identical results), and is now used by the tutorials. A tolerance of `1e-9` is normally sufficient.
* `nAlphaSubCycles`: reducing it from 3 to 1 made PowderBed2D 15 % faster but visibly changed the melt pool within 50 time steps, so it should only be changed after checking the results.
* `PIMPLE`: `nCorrectors 3` is usually enough.
* `MELTING`: the number of liquid-fraction corrector iterations is reported in the log each time step. If it regularly reaches `maxTempCorrector`, the liquid fraction is not converged to `epsilonTolerance`.
* Output: use `writeFormat binary;`, write less often and set `writeDiagnostics false;`.
* Parallel: the solver scales well; PowderBed3D (216k cells) ran at 1.88, 0.95 and 0.51 s per time step on 1, 2 and 4 cores. A decomposition with no split in the y (beam) direction keeps each beam column on one processor, which avoids communication in the heat source ray-tracing.

The heat source ray tracing walks precomputed, y-ordered columns of cells instead of searching the mesh for every cell each time step. Together with the removal of repeated work from the liquid-fraction corrector and PISO loops, this gave the following serial speed-ups over the previous version (OpenFOAM-13, identical results to within linear-solver round-off):

| Case | Cells | Previous (s/step) | Current (s/step) | Speed-up |
|---|---|---|---|---|
| PowderBed2D, melting (restart from t = 1e-5 s, 50 steps) | 115k | 48.3 | 1.58 | 30× |
| PowderBed2D, start-up (20 steps) | 115k | 46.8 | 1.0 | 47× |
| ArcCase (to t = 0.01 s, 1020 steps) | 7.2k | 0.272 | 0.066 | 4.1× |
| PowderBed3D (35 steps) | 216k | 2.84 | 2.21 | 1.3× |

The gain is largest for 2D (one-cell-thick) meshes, for which the mesh search used previously was particularly slow.

The implicit latent heat linearisation then reduced the number of liquid-fraction corrector iterations as follows (relative to the explicit update, with all of the above improvements in both):

| Case | Corrector iterations per step | Wall time | Result |
|---|---|---|---|
| PowderBed2D, melting (restart from t = 1e-5 s, 48 steps) | 26.8 → 5.1 | 80.3 s → 53.8 s (1.5×) | melt area within 0.01 % at t = 3e-5 s |
| GalliumCase (0–120 s) | 60.5 → 6.4 | 8050 s → 6969 s (1.2×) | melt fraction within 0.007 % and identical melt front at 60 s and 120 s |

With the explicit update, the GalliumCase liquid fraction frequently did not converge to `epsilonTolerance` within `maxTempCorrector` iterations; with the linearisation every time step converged in 7–19 iterations.

For the EB_3D laser welding case (864k cells, 4 cores, measured over short windows restarted from the same full-mesh state), `p_rgh` is now solved with `GAMG` instead of `PCG`/`DIC` (283 → about 20 iterations on the first corrector), with 3 instead of 5 PISO correctors (the 4th and 5th did no iterations), and `maxDeltaT` is raised from 2.5e-8 to 1e-7 s, so that the time step is limited by `maxCo 0.1` (to about 3e-8–8e-8 s):

| Window | Previous settings | Current settings | Speed-up | Result |
|---|---|---|---|---|
| Melting, t = 4e-5–4.5e-5 s | 1168 s | 243 s | 4.8× | peak temperature within 0.02 K, melt volume within 0.007 % |
| Vapourisation, t = 4.5e-5–7e-5 s | failed at t = 6.3e-5 s (5840 s estimated) | 1540 s | 3.8× | peak temperature and melt volume within 0.01 %, peak recoil pressure within 0.05 % of the linearised run at a time step of 2.5e-8 s |

The pressure solver and corrector changes on their own gave identical results (to 1e-7 relative) at 1.8× the speed. The case also now writes in binary. Its `endTime` of 2 s lies far beyond the t = 0.12 s at which the beam leaves the domain; reducing it is the largest remaining saving if only the weld itself is of interest.

`EB_3D_coarse` is a reduced version of EB_3D for workstations and laptops, run with `./Allrun parallel` on 12 processes. It uses the same set-up on a mesh coarsened by two in each direction (108k instead of 864k cells), with `maxCo 0.2`, `maxDeltaT 4e-7` and `endTime 0.12` s (when the beam leaves the domain). On 4 cores it took about 800 s per millisecond of simulated time once the keyhole had formed, so the whole weld takes about a day on 4 cores and correspondingly less on more. It was run stably to t = 10.5 ms. Compared with EB_3D:

* Heating is slower at first, since the beam energy is spread over larger cells (peak temperature 1928 K instead of 3195 K at t = 5e-5 s), but the peak temperature at vapourisation is the same (3749 K instead of 3752 K).
* The finer details of the melt pool and keyhole are lost, so it is suited to trying out the case and the parameters, not to quantitative results.

From t = 5 to 5.5 ms, with the original time step (`maxCo 0.1`, `maxDeltaT 1e-7`) as the reference, the current settings are 3.6× faster, with the melt volume within 0.2 %, the mean metal temperature within 0.01 K and the same peak temperature. The keyhole flow is chaotic, so individual cell values differ between any two time steps, including 1e-7 and 2e-7 s.

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


