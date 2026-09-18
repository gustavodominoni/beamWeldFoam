---
name: openfoam-solver-port
description: Port a custom OpenFOAM solver across major versions, especially from a standalone application solver (OpenFOAM 6-10 style, a main() built from #include headers) to a solver module run by foamRun (OpenFOAM 11+). Use this whenever someone wants to upgrade, port, modernise or migrate an OpenFOAM solver to a newer version, when a custom solver no longer compiles after an OpenFOAM upgrade, when converting case files between OpenFOAM versions (transportProperties to physicalProperties, turbulenceProperties to momentumTransport, cAlpha to interfaceCompression), or when they mention foamRun, solver modules, fvModels/fvConstraints, or an OpenFOAM version number alongside a solver name. Also use it for the reverse question, "what would it take to port X", since the survey steps answer that.
---

# Porting an OpenFOAM solver across major versions

Most of the work in a version port is *deletion*, not translation. A typical custom
solver is 20% original physics and 80% copied library internals that have drifted out
of date. The goal is to find the 20%, express it against the new framework, and throw
the rest away. A successful port usually has a strongly negative line count.

## The architectural shift in OpenFOAM 11+

Before OpenFOAM 11, a solver was an application: a `main()` assembled from textual
`#include` fragments (`createFields.H`, `UEqn.H`, `pEqn.H`) with the time loop written
out in full.

From OpenFOAM 11 the time loop lives in `foamRun`, and a solver is a *class* deriving
from a base module, overriding virtual functions that `foamRun` calls in order:

```
preSolve → moveMesh → motionCorrector → fvModels.correct → prePredictor
    → momentumTransportPredictor → momentumPredictor
    → thermophysicalPredictor → pressureCorrector
    → momentumTransportCorrector → postSolve
```

This matters because phase transport, PIMPLE control, Courant limiting, mesh motion and
turbulence all move into the base class. Anything the old solver carried locally to do
those jobs is now dead weight.

## Step 1: Establish ground truth before writing code

Never port from memory of the API. Methods get renamed, and some get deleted outright
with no replacement — if you guess, you will write plausible code that does not compile,
or worse, quietly restore something the framework removed on purpose.

Clone the target version's source. You need it to build anyway, and it is the only
authoritative API reference:

```bash
git clone --depth 1 --branch master https://github.com/OpenFOAM/OpenFOAM-13.git
```

Then do three readings before touching the solver:

1. **Pick the base module** and read it end to end. It is your template — you are
   writing a sibling of it. Map the old solver to its modern equivalent:

   | Old application | Base module |
   |---|---|
   | `interFoam`, `interDyMFoam` | `incompressibleVoF` |
   | `compressibleInterFoam` | `compressibleVoF` |
   | `multiphaseInterFoam` | `incompressibleMultiphaseVoF` |
   | `pimpleFoam` | `incompressibleFluid` |
   | `buoyantPimpleFoam`, `rhoPimpleFoam` | `fluid` / `isothermalFluid` |
   | `reactingFoam` | `multicomponentFluid` |
   | `laplacianFoam` on a solid | `solid` |

   Run `ls applications/modules/` for the full list. Read the base class chain too
   (`incompressibleVoF` → `twoPhaseVoFSolver` → `twoPhaseSolver` → `VoFSolver` →
   `fluidSolver` → `solver`) so you know what you inherit and what is protected.

2. **Read a tutorial of that module** under `tutorials/<module>/`. This is the case
   file format you must convert to, and it is faster than inferring it from code.

3. **Diff the old solver against the library version it was forked from.** Everything
   that matches upstream is deletable. Everything that differs is the physics you must
   preserve. This single step usually defines the whole job.

## Step 2: Map the old loop onto overrides

Work out where each old fragment lands before writing anything:

| Old fragment | New home |
|---|---|
| `createFields.H` | constructor initialiser list |
| `readControls.H`, custom dictionary reads | `read()` override |
| `alphaEqnSubCycle.H`, `alphaEqn.H` | delete — base class `prePredictor()` |
| property updates, source terms | `prePredictor()`, after calling the base |
| `UEqn.H` | `momentumPredictor()` |
| `TEqn.H` / energy / species | `thermophysicalPredictor()` |
| `pEqn.H` | `pressureCorrector()` |
| post-mesh-motion fixes | `motionCorrector()` |
| `CourantNo.H`, `setDeltaT.H`, `setRDeltaT.H` | delete — base class |
| `correctPhi.H`, `initCorrectPhi.H` | delete — base class |

Split the implementation across files by override rather than keeping one large `.C`.
It maps cleanly onto how the framework calls you and keeps each file reviewable.

## Step 3: Delete the vendored library code

Be aggressive here. Headers like `alphaEqn.H`, `alphaEqnSubCycle.H`, `alphaCourantNo.H`,
`createAlphaFluxes.H`, `setDeltaT.H`, `setRDeltaT.H`, `correctPhi.H`, `initCorrectPhi.H`,
`rhofs.H` and `alphaSuSp.H` are almost always verbatim copies of library internals. They
existed because the old architecture had no way to inherit them. Delete them and let the
base class do the work.

Keep a fragment only when you can point at the specific line where it differs from
upstream and explain why that difference is the solver's physics.

## Step 4: Translate the physics

Now port the code that is genuinely yours. Read `references/api-changes.md` for the
rename tables (`fvOptions` → `fvModels`/`fvConstraints`, `setRefCell` →
`pressureReference`, `polyMesh::findCell` → `meshSearch`, and the rest), and grep the
source tree whenever a symbol is not in the table.

Two habits prevent most of the pain:

- **Syntax-check continuously instead of building.** A full OpenFOAM build takes an hour;
  a syntax-only compile takes seconds and catches essentially every API error. Use
  `scripts/of-syntax-check.sh`. Only do a real build when you need to link and run.
- **Audit dimensions field by field against the original.** Write the old `dimensionSet`
  and the new named constant side by side and confirm they agree. Dimension errors are
  the easiest way to silently change physics, and OpenFOAM will happily run with a
  consistent but wrong set.

Preserve behaviour by default. A port is not the moment to improve the model — if you
spot a genuine bug, fix it, but list it explicitly so the owner can distinguish "ported"
from "changed".

## Step 5: Convert the case files

Script this rather than hand-editing. Cases share structure, hand edits drift, and a
script is re-runnable when you discover a mistake on case 9 of 10. See
`references/case-files.md` for the conversions and a converter to adapt.

Inventory any dictionary entry the new version no longer supports. Some removals change
results silently — these are the most important thing to tell the owner about, because
nothing will error and the numbers will simply differ.

## Step 6: Verify

Read `references/verification.md` for the full playbook. The short version, in
increasing order of strength:

1. Compiles clean with warnings on.
2. Every tutorial runs; dictionaries parse for any case too large to run.
3. Source terms cross-checked against their analytic form.
4. **Results compared against the original solver on the same case.**

Only the fourth actually demonstrates the port preserved the physics. The first three
are necessary and much cheaper, so do them first, but do not let them stand in for the
fourth. If you cannot do the comparison, say so plainly rather than implying the port is
validated.

## Step 7: Report what you did and did not establish

State which cases ran, which did not and why, which behaviours changed, and which
verification you were unable to perform. A port that quietly changes results is worse
than one that fails loudly, so the caveats are part of the deliverable, not an apology.

## Traps worth knowing in advance

- **A missing method may be deleted, not moved.** `MRF.correctBoundaryVelocity` no longer
  exists in OpenFOAM 13. Before "restoring" something, grep the whole source tree; if it
  is absent everywhere, the framework dropped it and your code should too.
- **`LIB_SRC` is empty** outside a `wmake` invocation. Use `$FOAM_SRC` in your own scripts.
- **`fvCFD.H` is not for modules.** Include what you use: `fvc.H`, `fvm.H`,
  `surfaceInterpolate.H`, `zeroGradientFvPatchFields.H`.
- **Constructor initialiser order must match declaration order**, and these classes have
  many members. Compile with `-Wall -Wextra` so reordering is caught.
- **Cached mesh geometry goes stale** if the mesh can move. Anything derived from
  `mesh.C()` or cell dimensions and stored in the constructor needs rebuilding in
  `motionCorrector()`, or the solver needs to reject moving meshes.
- **`read()` is re-invoked at runtime** through the object registry when `fvSolution`
  changes, so `runTimeModifiable` keeps working if you put dictionary reads there.
