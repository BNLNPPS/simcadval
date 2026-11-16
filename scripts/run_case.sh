#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# Full chain for one case, run inside the eic container (see ../eic.sh):
#   DD4hep compact --dd4hep2step--> STEP (+ materials.xml)
#   simcadval_scan --world dd4hep (reference: DDG4 Geant4Converter)
#   simcadval_scan --world step   (candidate: G4OCCTAssemblyVolume)
#   compare.py                    (metrics.json, plots, mismatch CSVs)
#
# Usage: run_case.sh NAME COMPACT [scan args...]
#   env STEP=<file>      use an existing STEP (+ STEP.materials.xml) instead of exporting
#   env EXPORT_ARGS=...  extra dd4hep2step args (e.g. --perturb name:dx,dy,dz)
#   env MATMAP=<file>    override the material map used by the step world
#   env SCANNER=<exe>    simcadval_scan binary (e.g. one linked against another G4OCCT)
#   env BUILD=<dir>      build directory with dd4hep2step and simcadval_scan (default ./build)
#   env OUT=<dir>        output root (default ./out); results go to OUT/NAME
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
name=$1 compact=$2
shift 2
out=${OUT:-$here/out}/$name
mkdir -p "$out"
bin=${BUILD:-$here/build}

step=${STEP:-$out/$name.step}
if [ -z "${STEP:-}" ]; then
  # shellcheck disable=SC2086
  "$bin/dd4hep2step" "$compact" "$step" --report "$out/export.json" ${EXPORT_ARGS:-} 2>&1 \
    | grep -E "^dd4hep2step|volume check|FAILED|read-back" | tee "$out/export.log"
fi
report=$out/export.json
[ -f "$report" ] || report=$(dirname "$step")/export.json

for w in dd4hep step; do
  extra=()
  [ "$w" = step ] && extra=(--step "$step" ${MATMAP:+--matmap "$MATMAP"})
  rm -f "$out/$w.hung" "$out/scan_$w.log"
  first=0
  while :; do
    rc=0
    "${SCANNER:-$bin/simcadval_scan}" --compact "$compact" --world $w "${extra[@]}" --out "$out/$w" "$@" \
      --first-ray $first >>"$out/scan_$w.log" 2>&1 || rc=$?
    [ $rc -eq 0 ] && break
    if [ $rc -eq 3 ] && tail -1 "$out/$w.hung" | grep -q '^ray '; then
      first=$(( $(tail -1 "$out/$w.hung" | awk '{print $2}') + 1 ))   # resume after the hung ray
      echo "  [$w] ray $((first - 1)) hung; resuming at $first"
      continue
    fi
    echo "scan $w failed (rc=$rc), see $out/scan_$w.log"; tail -5 "$out/scan_$w.log"; exit 1
  done
  grep -E "^simcadval_scan|G4Exception |killed" "$out/scan_$w.log"
done
python3 "$here/scripts/compare.py" "$out/dd4hep" "$out/step" "$out/cmp" --label-a dd4hep --label-b step \
  --export-report "$report" | tee "$out/summary.txt"
