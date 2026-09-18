#!/bin/bash
#------------------------------------------------------------------------------
# Syntax-only compile of an OpenFOAM solver module against the target version's
# headers. A full OpenFOAM build takes about an hour; this takes seconds and
# catches essentially every API error, so it is the right inner loop while
# porting. Run a real build only when you need to link and execute.
#
# Usage:
#   of-syntax-check.sh <solver-dir> [file.C ...]
#
# Expects the OpenFOAM environment to be sourced already. Note that LIB_SRC is
# only set by wmake itself, so we fall back to FOAM_SRC.
#------------------------------------------------------------------------------
set -u

SOLVER_DIR=${1:?usage: of-syntax-check.sh <solver-dir> [file.C ...]}
shift || true

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    echo "Source the OpenFOAM etc/bashrc first." >&2
    exit 1
fi

LIB_SRC=${LIB_SRC:-$FOAM_SRC}

# Link headers for any library not yet built, so includes resolve.
for d in applications/modules/* src/*; do
    [ -d "$WM_PROJECT_DIR/$d/Make" ] && [ ! -d "$WM_PROJECT_DIR/$d/lnInclude" ] \
        && wmakeLnInclude "$WM_PROJECT_DIR/$d" > /dev/null 2>&1
done

INC=""
for d in $(ls -d "$FOAM_MODULES"/*/lnInclude 2>/dev/null) \
         $(ls -d "$LIB_SRC"/*/lnInclude "$LIB_SRC"/*/*/lnInclude 2>/dev/null); do
    INC="$INC -I$d"
done
INC="$INC -IlnInclude -I. -I$LIB_SRC/OpenFOAM/lnInclude -I$LIB_SRC/OSspecific/POSIX/lnInclude"

cd "$SOLVER_DIR" || exit 1
FILES=${*:-$(ls ./*.C 2>/dev/null)}

rc=0
for f in $FILES; do
    if g++ -std=c++14 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 \
           -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter \
           -Wno-invalid-offsetof -Wno-attributes -O0 -DNoRepository \
           -ftemplate-depth-256 $INC -fPIC -fsyntax-only "$f"
    then
        echo "ok    $f"
    else
        echo "FAIL  $f"
        rc=1
    fi
done
exit $rc
