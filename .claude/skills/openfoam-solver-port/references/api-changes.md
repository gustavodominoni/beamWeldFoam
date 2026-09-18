# API changes: OpenFOAM 6 → 13

Verified against the OpenFOAM-13 source during a real `interFoam`-derived port. Where a
symbol is not listed, grep the target source tree rather than guessing — this table is a
head start, not a complete catalogue.

## Contents

- [Solver structure](#solver-structure)
- [Sources and constraints](#sources-and-constraints)
- [Fields and mesh access](#fields-and-mesh-access)
- [Dictionaries and lookups](#dictionaries-and-lookups)
- [Pressure, gravity and reference cells](#pressure-gravity-and-reference-cells)
- [Turbulence](#turbulence)
- [Removed with no replacement](#removed-with-no-replacement)
- [Case dictionary renames](#case-dictionary-renames)

## Solver structure

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `int main()` + `#include` fragments | class deriving from a module in `applications/modules/` |
| `#include "fvCFD.H"` | `fvc.H`, `fvm.H`, `surfaceInterpolate.H` as needed |
| `Make/files` → `EXE = $(FOAM_USER_APPBIN)/x` | `LIB = $(FOAM_USER_LIBBIN)/libx` |
| `while (runTime.run())` written out | `foamRun` owns the loop; you override hooks |
| `pimple.loop()` in your code | base class |
| registration: none | `TypeName("x")` + `addToRunTimeSelectionTable(solver, x, fvMesh)` |

A module is selected by `solver x;` in `controlDict`, or `foamRun -solver x`. `foamRun`
loads `libx.so` by name, so the library name and the registered type name must match.

## Sources and constraints

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `fvOptions(rho, U)` | `fvModels().source(rho, U)` |
| `fvOptions.constrain(eqn)` | `fvConstraints().constrain(eqn)` |
| `fvOptions.correct(U)` | `fvConstraints().constrain(U)` |
| `#include "createFvOptions.H"` | inherited; call `fvModels()` / `fvConstraints()` |

## Fields and mesh access

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `runTime.timeName()` | `runTime.name()` |
| `mesh.solutionDict()` | `mesh.solution()` |
| `mesh.ddtScheme("ddt(alpha)")` | `mesh.schemes().ddt("ddt(alpha)")` |
| `mesh.setFluxRequired(p.name())` | `mesh.schemes().setFluxRequired(p.name())` |
| `field.internalField()` (non-const) | `field.primitiveFieldRef()` |
| `field.ref()` on the internal field | `field.internalFieldRef()` |
| `mesh.findCell(pt)` | `meshSearch::New(mesh).findCell(pt)` |
| `vector::zero` | `Zero` |

`meshSearch::New(mesh)` returns a cached octree registered on the mesh, so it is cheaper
than the old linear search and is invalidated automatically on mesh change. Include
`meshSearch.H`. It returns `-1` for a point outside the mesh — always guard, since
indexing a field with `-1` is undefined behaviour.

## Dictionaries and lookups

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `readScalar(dict.lookup("k"))` | `dict.lookup<scalar>("k")` |
| `readLabel(dict.lookup("n"))` | `dict.lookup<label>("n")` |
| `dimensionedScalar("k", dims, dict.lookup("k"))` | `dimensionedScalar("k", dims, dict.lookup<scalar>("k"))` |
| `dict.lookupOrDefault<T>("k", v)` | unchanged |

Named dimension sets read better than raw exponents and are checked by the compiler:
`dimSpecificHeatCapacity`, `dimThermalConductivity`, `dimEnergy/dimMass`, `dimPressure`,
`dimDensity`, `dimTemperature`, `dimPower/dimVolume`, `dimMass/dimMoles`.

## Pressure, gravity and reference cells

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `setRefCell(p, p_rgh, dict, cell, value)` | `pressureReference pr(p, p_rgh, dict)` |
| `getRefCellValue(p, pRefCell)` | `getRefCellValue(p, pr.refCell())` |
| `pRefValue` | `pr.refValue()` |
| `#include "readGravitationalAcceleration.H"` | `solvers::buoyancy` member of `VoFSolver` |
| `gh`, `ghf`, `ghRef` | `buoyancy.gh`, `buoyancy.ghf`, `buoyancy.ghRef` |

The pressure equation also changed sign convention. OpenFOAM 6:

```cpp
fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
phi = phiHbyA - p_rghEqn.flux();
```

OpenFOAM 13:

```cpp
fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh) == source
phi = phiHbyA + p_rghEqn.flux();
```

Copy the convention from the base module wholesale. Mixing the two produces a solver
that converges to the wrong pressure field, which is easy to miss.

## Turbulence

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `turbulence->divDevRhoReff(rho, U)` | `divDevTau(U)` (virtual on the module) |
| `turbulence->correct()` | `momentumTransportCorrector()` override |
| `turbulence->validate()` | handled by the base class |
| `incompressible::turbulenceModel` | `incompressibleInterPhaseTransportModel` for VoF |

## Removed with no replacement

Check for these before assuming an omission is a mistake:

- `MRF.correctBoundaryVelocity(U)` — gone entirely. MRF boundary handling now goes
  through `makeRelative`/`makeAbsolute` and `MRF.update()`, which the base classes call.
- `cAlpha`, `scAlpha`, `icAlpha` in `fvSolution` — replaced by run-time selectable
  schemes on `div(phi,alpha)` in `fvSchemes`. `cAlpha` maps onto the scheme's coefficient,
  but **`scAlpha` (shear compression) and `icAlpha` (isotropic compression) have no
  equivalent**. Cases that used them will produce different interface sharpening. Flag
  this to the owner; it is a real behavioural change, not a cosmetic one.
- `Time::timeName()` as an instance method — use `name()`.

`twoPhaseVoFSolver` raises a fatal error if `cAlpha` is still present in the alpha solver
dictionary, so leftover entries fail loudly rather than being ignored.

## Case dictionary renames

| OpenFOAM 6 | OpenFOAM 13 |
|---|---|
| `constant/transportProperties` | `constant/phaseProperties` + `constant/physicalProperties.<phase>` |
| `constant/turbulenceProperties` | `constant/momentumTransport` |
| `constant/motionProperties` | `constant/dynamicMeshDict` (or omit for a static mesh) |
| `application foo;` in `controlDict` | `solver foo;` |
| `phases (a b);` + `a { ... }` sub-dicts | `phases (a b);` in `phaseProperties`, per-phase files |
| `transportModel Newtonian;` | `viscosityModel constant;` |

Per-phase files hold `viscosityModel`, `nu`, `rho` plus any custom properties. The mixture
file holds `phases`, `sigma` and custom mixture properties. A custom solver reads its own
entries from whichever of these it likes — `viscosityModel` derives from `IOdictionary`,
so `mixture.nuModel1().lookup<scalar>("cp")` works directly.
