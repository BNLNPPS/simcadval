---
name: investigate-mismatch
description: Diagnose why the DD4hep and STEP worlds of a simcadval case disagree. Use when asked to investigate, explain or diagnose a mismatch in an out/NAME case (after run_case.sh), or what triage found. Replays rays, runs targeted rescans and writes out/NAME/cmp/diagnosis.md with the evidence.
---

# Investigate a simcadval mismatch

Goal: for one case directory `out/NAME` (made by `scripts/run_case.sh`), say
**what** differs between World A (DD4hep) and World B (STEP via G4OCCT),
**where** (which parts), **by how much**, and show the evidence. Write the
result to `out/NAME/cmp/diagnosis.md`.

Scope: diagnose the comparison; do not debug G4OCCT. When the cause is G4OCCT
navigation, classify it as a navigation artifact, give the evidence and stop.
Do not edit sources, examples or G4OCCT, and do not build other G4OCCT
versions. Write only under `out/NAME/` (rescans and scratch scripts go in
`out/NAME/investigation/`).

Run everything from the repo root inside the container: wrap every command,
including the ones triage prints bare (`python3 …`, `${SCANNER:-…} …`), as
`./eic.sh '<command>'` with single quotes. If `build/` has no
`simcadval_scan`, build first (README, Quick start).

## Categories

Every finding gets exactly one:

| Category | What it means | Typical evidence |
|---|---|---|
| placement | a part is shifted, rotated or the wrong size in the CAD | same material sequence, one boundary offset at that part by an amount consistent across rays (see step 7); with no other family present also a low IoU for that part and boundary F1 near 1 |
| material | a part has the wrong material in World B | identical step boundaries; the part's material differs on whole steps; points disagree with the same pair (e.g. Kapton -> Copper); volume-id agreement unaffected; the map entry for that part differs from its DD4hep material |
| missing part | a part is absent from World B | the part's material never appears in B along rays or at points where A has it; `--debug-ray` finds no part there; check `export.json` status |
| navigation artifact | the STEP shapes are right but G4OCCT navigates them wrongly | rays disagree while points agree; whole steps lost on entry through a face (often B reads the world material, e.g. Air); step boundaries unchanged. Known case: thin touching slabs, eic/G4OCCT#405 |
| no difference | worlds agree at floating-point level | 1 - material agreement below about 1e-8, no point disagreements, boundary F1 @1nm = 1 |
| unclear | evidence does not separate the categories | say what was tried and what would settle it |

One case can hold several families, e.g. a CAD fault on top of the navigation
artifact. Diagnose each family separately. Global numbers (worst volumes,
boundary F1, max boundary shift, nearest-boundary percentiles) mix all
families, so a large family hides a small one there: measure each family on
its own rays.

## Procedure

1. **Read the headline and the case's files.** `out/NAME/summary.txt`, then
   `out/NAME/cmp/metrics.json`: material agreement, point disagreements and
   their pairs, worst volumes, and robustness fields (`rays_missing`, `hung`,
   `killed_events`, `g4_exceptions` under `a`/`b`). Check `export.log` and
   `export.json`: a `perturb` entry means the STEP was moved on purpose by
   `dd4hep2step --perturb`; report it. Other comparison outputs in the case
   directory (e.g. `cmp_*`) are valid evidence; name them if you use them.
   If the "no difference" row fits, write that and stop.
2. **Triage.** If `out/NAME/cmp/triage/triage.md` is missing, run
   `python3 scripts/triage.py out/NAME`. Its strata table lists every family
   (reference material -> candidate material, steps, share of length),
   including tiny ones. Its failure modes cover only the large families; each
   has a representative ray with ready `show_ray.py` and `--debug-ray`
   commands.
3. **Get a ray for every family**, including strata without a mode. For
   those, take rays from `out/NAME/cmp/triage/features.csv.gz` (one row per
   sampled step): rows with `label` = 1 (mismatching) and the stratum's
   `material` -> `material_b`; `evt` is the ray, `s0_mm` and `len_mm` the
   step, `volume`, `r`, `z`, `phi` where it is. (`cmp/ray_mismatches.csv`
   keeps only the 5000 longest intervals, so small families are usually
   missing from it.) Replay with
   `python3 scripts/show_ray.py out/NAME/dd4hep out/NAME/step EVT SMIN SMAX --report out/NAME/export.json`
   (s range = the step ± 1 mm). Note whether boundaries match and only the
   material differs, whether a boundary is offset (by how much), or whether a
   part is absent.
4. **Ask the parts.** Run the `--debug-ray` command (from triage, or the same
   form with your ray and s range). It writes no files. Geant4 prints about
   1400 lines first; keep the useful part with `| awk '/^s=/{p=1} p'`. It
   lists the STEP parts near each step point and what their solids answer
   (`Inside`, safety, distances in and out). A part answering `outside` with
   safety 0 on its own face also happens where B enters correctly, so compare
   a failing ray with one that enters the same kind of part correctly before
   drawing a conclusion.
5. **Check the material map** for any material finding. The map is
   `<STEP>.materials.xml` (STEP path in `out/NAME/step.meta.json`), or
   `out/NAME/materials.xml` when one exists there (the case was run with
   `MATMAP=`). Compare the part's entry (`stepName="<part>"`) with its
   DD4hep material in `out/NAME/export.json` (or the export report next to
   the STEP file).
6. **Confirm with a targeted rescan** when one ray is not enough (to measure
   a shift, or to check that a pattern holds across a part). Aim rays at the
   part, reusing the case's STEP file and map:
   `OUT=out/NAME/investigation STEP=<step> [MATMAP=<map>] scripts/run_case.sh zoomK <compact> --rays 10000 --ray-mode aim --box <hx,hy,hz> --target <cx,cy,cz,hx,hy,hz> --points 0`.
   In aim mode `--box` (half-lengths, mm) is where rays start and `--target`
   (centre and half-lengths, mm) is the box they aim at; take the target from
   replay positions or the part's placement. Keep each rescan at 20 000 rays
   or fewer. With `--points 0`, compare.py prints "Mean of empty slice"
   warnings and `nan` point agreement; ignore them.
7. **Measure magnitudes per family.** For placement: a face with unit normal
   n moved by a displacement D shifts the boundary along a ray with direction
   u by Δs = (D·n)/(u·n), so Δs·(u·n) is the same for every ray through that
   face. Compute it over many rays per face (from a rescan) rather than from
   the global max boundary shift, which a lost boundary or a corner clip
   inflates. Small analysis scripts may go in `out/NAME/investigation/`;
   read scans with `numpy.fromfile(prefix + ".steps.bin", dtype=STEP_DT)`
   and `RAY_DT`, both defined in `scripts/compare.py` (material names are
   in `<prefix>.meta.json`). For material: the affected length and the map
   entry. For a navigation artifact: the share of the mismatch length and
   the solid types and thicknesses involved.
8. **Write `out/NAME/cmp/diagnosis.md`:**
   - first line: the verdict in one sentence;
   - a table with one row per finding: category, parts (`<volume>@<id>`),
     magnitude, share of mismatch length, confidence (high / medium / low);
   - per finding, the evidence: each command run and the output lines that
     show it;
   - what is not explained, and what would settle it.

## Rules

- Every claim cites output from a command you ran in this investigation.
  Triage rules and modes are correlations: confirm them with a replay.
- Do not report a cause you have not shown. Use "unclear" instead.
- Keep G4OCCT internals out of the diagnosis beyond the observed answers of
  its solids; for a navigation artifact, point to a matching known issue
  (eic/G4OCCT#405) or say that it looks new.
