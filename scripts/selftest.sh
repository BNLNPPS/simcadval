#!/usr/bin/env bash
# Metric self-test: deliberately wrong CAD must be detected and localised.
# Run inside the eic container (see ../eic.sh) after building.
#
# Each faulty STEP world is compared with the STEP world of the correct export,
# both navigated by the same G4OCCT, so whatever G4OCCT itself gets wrong
# against DD4hep is identical in both and cancels; only the injected fault is left.
#   shift10um, shift1um  inner vertex-barrel modules moved 10 um / 1 um along z
#   matswap              one ECal Kapton slice prototype (slice3@18) read as Copper
# Null test: the DD4hep scans of two cases come from separate processes and must agree exactly.
#
# Usage: selftest.sh [scan args]   (default: 20000 ip rays, 10^6 points in the MiniSiD detector box)
#   env OUT=<dir> (default ./out; results in OUT/selftest), BUILD=<dir> (default ./build)
# Exit status 0 if every check passes, 1 otherwise.
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
root=${OUT:-$here/out}/selftest
bin=${BUILD:-$here/build}
C=$here/examples/MiniSiD.xml
if [ $# -gt 0 ]; then ARGS=("$@"); else
  ARGS=(--rays 20000 --ray-mode ip --points 1000000 --box 1450,1450,1800 --seed 12345); fi
export OUT=$root BUILD=$bin

mkdir -p "$root/base" "$root/matswap"
base=$root/base/base.step
"$bin/dd4hep2step" "$C" "$base" --report "$root/base/export.json" > "$root/base/export.log" 2>&1 || {
  echo "export failed, see $root/base/export.log"; exit 1; }
sed 's|stepName="slice3@18" geant4Name="Kapton"|stepName="slice3@18" geant4Name="Copper"|' \
  "$base.materials.xml" > "$root/matswap/materials.xml"
[ "$(grep -c 'stepName="slice3@18" geant4Name="Copper"' "$root/matswap/materials.xml")" = 1 ] || {
  echo "matswap: slice3@18 Kapton entry not found in $base.materials.xml"; exit 1; }

echo "selftest: scanning 4 cases in parallel (${ARGS[*]})"
pids=()
STEP=$base "$here/scripts/run_case.sh" base "$C" "${ARGS[@]}" > "$root/base.log" 2>&1 & pids+=($!)
EXPORT_ARGS="--perturb VtxBarrelModuleInner:0,0,0.010" "$here/scripts/run_case.sh" shift10um "$C" "${ARGS[@]}" \
  > "$root/shift10um.log" 2>&1 & pids+=($!)
EXPORT_ARGS="--perturb VtxBarrelModuleInner:0,0,0.001" "$here/scripts/run_case.sh" shift1um "$C" "${ARGS[@]}" \
  > "$root/shift1um.log" 2>&1 & pids+=($!)
STEP=$base MATMAP=$root/matswap/materials.xml "$here/scripts/run_case.sh" matswap "$C" "${ARGS[@]}" \
  > "$root/matswap.log" 2>&1 & pids+=($!)
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ $rc -eq 0 ] || { echo "a scan failed, see $root/*.log"; exit 1; }

cmp_case() {  # cmp_case NAME A_PREFIX B_PREFIX LABEL_B
  python3 "$here/scripts/compare.py" "$2" "$3" "$root/$1/cmp_$4" --label-a "${2##*/}" --label-b "$4" \
    --export-report "$root/base/export.json" --no-plots > "$root/$1/cmp_$4.txt"
}
for f in shift10um shift1um matswap; do cmp_case "$f" "$root/base/step" "$root/$f/step" vs_base; done
cmp_case shift10um "$root/base/dd4hep" "$root/shift10um/dd4hep" null

python3 - "$root" <<'EOF'
import json, sys
root = sys.argv[1]
def load(case, tag):
    return json.load(open(f"{root}/{case}/cmp_{tag}/metrics.json"))
fails = []
def check(name, ok, detail):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}")
    if not ok:
        fails.append(name)

print("null test (DD4hep world, two processes)")
m = load("shift10um", "null")
check("null: identical", m["rays"]["material_mismatch_length_mm"] == 0 and m["headline"]["max_boundary_shift_mm"] == 0
      and m["points"]["n_disagree"] == 0,
      f"mismatch {m['rays']['material_mismatch_length_mm']} mm, max shift {m['headline']['max_boundary_shift_mm']} mm")

mis = {}
for case, dz in (("shift10um", 0.010), ("shift1um", 0.001)):
    print(f"{case} (vs correct STEP)")
    m = load(case, "vs_base")
    mis[case] = m["rays"]["material_mismatch_length_mm"]
    w = m["worst_volumes"]
    moved = {"component0", "component1"}
    in_moved = sum(v["mismatch_mm"] for v in w if v["name"] in moved)
    elsewhere = sum(v["mismatch_mm"] for v in w if v["name"] not in moved)
    check(f"{case}: detected", mis[case] > 0, f"material mismatch {mis[case]:.4g} mm")
    check(f"{case}: localised", {v["name"] for v in w[:2]} == moved and elsewhere <= 1e-6 * in_moved,
          f"top volumes {[v['name'] + '@' + str(v['id']) for v in w[:2]]}, mismatch elsewhere {elsewhere:.3g} mm")
    s = m["headline"]["max_boundary_shift_mm"]
    check(f"{case}: size", 0.5 * dz <= s <= 2 * dz, f"max boundary shift {s:.4g} mm for an injected {dz} mm")
ratio = mis["shift10um"] / mis["shift1um"] if mis["shift1um"] > 0 else float("inf")
check("shift: scales with the fault", 5 <= ratio <= 20, f"mismatch ratio 10um/1um = {ratio:.2f}")

print("matswap (vs correct STEP)")
m = load("matswap", "vs_base")
bad = sorted(k for k, v in m["rays"]["material_iou"].items() if v["iou"] < 1 - 1e-12)
check("matswap: detected", m["rays"]["material_mismatch_length_mm"] > 0,
      f"material mismatch {m['rays']['material_mismatch_length_mm']:.4g} mm")
check("matswap: geometry untouched", m["rays"]["volume_mismatch_fraction_nonworld"] == 0,
      f"volume mismatch {m['rays']['volume_mismatch_fraction_nonworld']}")
check("matswap: only Kapton/Copper", bad == ["Copper", "Kapton"], f"materials with IoU < 1: {bad}")
p = m["points"]
if p["n"]:
    check("matswap: points", p["n_disagree"] > 0 and set(p["disagreements"]) == {"Kapton -> Copper"},
          f"{p['n_disagree']} points disagree: {p['disagreements']}")

print(f"selftest: {'all checks passed' if not fails else 'FAILED: ' + ', '.join(fails)}")
sys.exit(1 if fails else 0)
EOF
