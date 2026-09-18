# Converting case files between OpenFOAM versions

Write a converter script rather than editing cases by hand. Cases in a tutorial suite
share structure, hand edits drift apart, and you will almost always discover a mistake
while converting the last case that has to be applied to all the earlier ones.

## What the converter has to do

Working from an OpenFOAM 6 VoF case to OpenFOAM 13:

1. **Split `constant/transportProperties`.** The `phases (a b)` list plus the per-phase
   sub-dictionaries become `constant/phaseProperties` (phase list, `sigma`, mixture-level
   custom entries) and one `constant/physicalProperties.<phase>` per phase (`viscosityModel`,
   `nu`, `rho`, per-phase custom entries). `transportModel Newtonian;` becomes
   `viscosityModel constant;`.
2. **Rename `constant/turbulenceProperties`** to `constant/momentumTransport`, carrying
   `simulationType` across.
3. **Delete dictionaries the new solver never reads** — `motionProperties`, stale
   `thermophysicalProperties.*` left over from another solver.
4. **Swap `application foo;` for `solver foo;`** in `controlDict`.
5. **Move interface compression** from `cAlpha` in `fvSolution` to the `div(phi,alpha)`
   scheme in `fvSchemes`:
   `div(phi,alpha) Gauss vanLeer;` + `cAlpha 1;` becomes
   `div(phi,alpha) Gauss interfaceCompression vanLeer 1;`
   Delete `div(phirb,alpha)`, which no longer exists.
6. **Normalise types.** `pRefCell 0.0;` must become `pRefCell 0;` — it is read as a
   `label` and a decimal point is a parse error.
7. **Add `Allrun`/`Allclean`** if absent, calling `foamRun -solver <name>`.

## Approach that works

Parse with a light regex pass rather than a full dictionary parser. Case dictionaries are
regular enough, and the alternative is a large dependency for a one-off job. Two habits
keep it safe:

- **Assert before replacing.** Check each pattern matches exactly once before substituting,
  so a case with unexpected formatting fails loudly instead of being silently skipped.
- **Capture values before deleting them.** `cAlpha` has to be read out of `fvSolution`
  before that line is removed, because it becomes the scheme coefficient in `fvSchemes`.

```python
# Read cAlpha first, then delete the compression controls, then use the captured
# value when rewriting the div(phi,alpha) scheme.
m = re.search(r"^\s*cAlpha\s+([0-9.eE+-]+)\s*;", text, flags=re.M)
cAlpha = m.group(1) if m else "1"
text, n = re.subn(r"^[ \t]*(cAlpha|scAlpha|icAlpha)\s+[^;]+;[^\n]*\n", "", text, flags=re.M)

text, n = re.subn(
    r"^([ \t]*div\(phi,alpha\))\s+Gauss\s+(\w+)\s*;[^\n]*$",
    rf"\1  Gauss interfaceCompression \2 {cAlpha};",
    text, flags=re.M)
assert n == 1, "div(phi,alpha) scheme not found"
```

## Inventory what you drop

Keep a record of every entry the new version cannot express, and report it. Most removals
are harmless, but some change results with no error — `scAlpha` and `icAlpha` are the
known examples for VoF. The owner needs to know which of their cases were affected,
because nothing downstream will warn them.

Counting removals per case as the script runs is enough to catch this: a case that loses
two compression controls instead of one had `scAlpha` set.

## Validate the conversion

Parsing is a cheap, complete check that scales to cases too large to run:

```bash
for f in constant/* system/*; do foamDictionary "$f" > /dev/null || echo "FAIL $f"; done
```

Run it across every case and every `0/` field file. It catches truncation, unbalanced
braces and type errors in seconds, and it is the only verification available for a case
whose mesh will not fit in memory.

`foamDictionary` is not built by a partial OpenFOAM build — build it explicitly with
`wmake applications/utilities/miscellaneous/foamDictionary` if it is missing, and do not
mistake "command not found" for a parse failure.
