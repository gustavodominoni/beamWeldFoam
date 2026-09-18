# Verifying a solver port

The question a port has to answer is "does this compute the same physics as before".
Compilation and clean runs do not answer it. They are cheap and necessary, so do them
first, but treat them as prerequisites rather than evidence.

## The ladder, weakest to strongest

### 1. Compiles clean

Build with `-Wall -Wextra` and read the warnings. In a port these are unusually
informative, because reordered constructor initialiser lists, unused variables left from
deleted code and sign-conversion surprises all show up here.

### 2. Cases run

Run every tutorial for a handful of time steps. OpenFOAM enables floating-point exception
trapping by default, so a run that reaches `End` has already ruled out division by zero
and NaN propagation through every code path it touched.

Check the run actually exercised the physics rather than trivially succeeding. A case
where the melt-fraction residual stays exactly zero, or the source field is uniformly
zero, has not tested the code you care about — find a case where the term is active.

### 3. Dictionaries parse

For cases too large to run, `foamDictionary` over every file at least proves the
conversion is syntactically valid. Weak, but it scales and it is better than nothing.

### 4. Source terms cross-checked analytically

This is the first step that says anything about correctness. For any term with a closed
form — a Gaussian heat source, an analytic boundary profile — compute the expected value
independently and compare against the field the solver wrote.

Two things make it sharp:

- **Check the discrete structure, not just the magnitude.** If a ray-tracing heat source
  should deposit on exactly one cell per column, count the non-zero cells and compare
  against the mesh dimensions. An exact integer match is strong evidence the indexing is
  right, and it is the kind of error a magnitude check alone will not catch.
- **Account for discretisation before calling a mismatch a bug.** A peak sampled at a cell
  centre offset from the true peak will read slightly low. Predict the offset and confirm
  it explains the residual, rather than dismissing the difference as rounding.

### 5. Compared against the original solver

The only step that demonstrates the port preserved the physics. Build the old solver
against its old OpenFOAM, run the same case in both, and compare the fields that matter.
For a validation case with published reference data, compare both against that data.

This is the step most likely to get skipped, because it needs two OpenFOAM installations
side by side. Skipping it is defensible; implying it was done is not. If you cannot run
it, say which cases remain unvalidated and what someone else would need to do.

## Dimensional audit

Before trusting any run, put every field and constant next to its original and confirm
the dimensions agree:

| Quantity | OpenFOAM 6 `dimensionSet` | Named equivalent |
|---|---|---|
| specific heat | `(0 2 -2 -1 0)` | `dimSpecificHeatCapacity` |
| conductivity | `(1 1 -3 -1 0)` | `dimThermalConductivity` |
| latent heat | `(0 2 -2 0 0)` | `dimEnergy/dimMass` |
| pressure | `(1 -1 -2 0 0)` | `dimPressure` |
| volumetric source | `(1 -1 -3 0 0)` | `dimPower/dimVolume` |
| surface flux | `(1 0 -3 0 0)` | `dimPower/dimArea` |
| molar mass | `(1 0 0 0 -1)` | `dimMass/dimMoles` |

A dimensionally consistent but wrong set runs happily and produces wrong numbers, so this
is worth doing by hand for every quantity rather than trusting the translation.

## Reporting

Separate what you established from what you assumed:

- which cases ran to completion, and which did not, with the actual reason
- resource failures distinguished from solver failures — an out-of-memory kill with no
  OpenFOAM error in the log is a machine limit; confirm it from the kernel log rather
  than inferring it
- behavioural changes, especially silent ones from dropped dictionary entries
- bugs fixed in passing, listed so the owner can tell "ported" from "changed"
- verification not performed

A port whose caveats are documented can be trusted incrementally. One presented as fully
validated when it is not will be found out on the first paper that uses it.
