---
name: investigate-mismatch
description: Diagnose why the DD4hep and STEP worlds of a simcadval case disagree. Use when asked to investigate, explain or diagnose a mismatch in an out/NAME case (after run_case.sh), or what triage found. Uses the svtool investigation tools and submits a checked diagnosis (out/NAME/cmp/diagnosis.json and diagnosis.md).
---

# Investigate a simcadval mismatch

Goal: for one case directory `out/NAME` (made by `scripts/run_case.sh`), say
**what** differs between World A (DD4hep, the reference) and World B (STEP
via G4OCCT), **where** (which parts), **by how much**, with evidence, and
submit it as a diagnosis.

## Tools

All work goes through the svtool tools. Each returns one JSON object and is
logged with a call id (`c001`, `c002`, ...); the diagnosis cites these ids.

- From a shell (Claude Code): `scripts/svtool.py TOOL --arg value ...` from the
  repo root (it enters the container itself). Lists: `--parts a,b`; vectors:
  `--target cx,cy,cz,hx,hy,hz`; the diagnosis: `--diagnosis @file.json`.
- As functions (agent harness): the same names and arguments;
  `scripts/svtool.py spec` prints the definitions.

| Tool | Use |
|---|---|
| `case_summary` | headline, worst parts, robustness, export status, `no_difference_test` |
| `mismatch_families` | every mismatch family with share, parts, region, examples and `boundary_pattern` |
| `point_mismatches` | points with different materials, grouped |
| `ray` | one ray, both worlds side by side, with which world has a boundary where |
| `debug_ray` | what each nearby STEP part answers along one ray |
| `material_map` | STEP material map vs DD4hep materials |
| `part_info` | a part's export record and placements (centres and sizes for rescans) |
| `boundary_shifts` | how far a part's boundaries moved; fits the displacement D |
| `rescan` | aimed rescan of one region; the result is a case for every tool |
| `triage` | statistical rules and failure modes (slow the first time) |
| `submit_diagnosis` | check and write the diagnosis |

Scope: diagnose the comparison, do not debug G4OCCT. Do not edit sources,
examples or G4OCCT. Everything you produce goes through the tools.

## Categories

Each finding gets exactly one category. The `boundary_pattern` of a
mismatch family (from `mismatch_families`) is the first sorting signal:

| Category | Pattern | Confirm with |
|---|---|---|
| `navigation_artifact`: B's shapes are right, G4OCCT navigates them wrongly | `same_boundaries`; B has the world material where A has a part's material; points agree | `debug_ray` on an example: the part answers `outside` with safety 0 on its own entry face; `material_map` entry correct. Known case: thin touching slabs, eic/G4OCCT#405 |
| `material`: a part has the wrong material in B | `same_boundaries`; B has another real material | `material_map` differences; `point_mismatches` show the same pair |
| `placement`: a part is shifted, rotated or the wrong size in B | `one_end_each` | `boundary_shifts` for the part: D from the fit |
| `missing_part`: a part of A is absent in B | `only_in_a` | `part_info` (export status, read-back); `debug_ray` finds no part there |
| `extra_part`: B has a part A lacks | `only_in_b` | `debug_ray` names the part |
| `no_difference` | `case_summary` `no_difference_test.passes` | nothing else needed |
| `unclear` | `mixed`, or evidence that does not separate categories | say what would settle it |

A case can hold several families with different causes (e.g. a CAD fault
on top of the navigation artifact). Diagnose each family on its own: the
global numbers in `case_summary` mix them, and a large family hides a small
one there.

## Procedure

1. `case_summary`. If `no_difference_test.passes`, submit one
   `no_difference` finding (parts `[]`, share 0) and stop. Note robustness
   problems (hung or killed rays, Geant4 exceptions) and failed export parts.
   A `perturb` entry in `export` means the STEP was moved on purpose; report
   it.
2. `mismatch_families` (default `min_mm` 1e-6 mm). Every family needs an
   explanation, including tiny ones. Families that mirror each other
   (`Carbon -> Air` in a part and `Air -> Carbon` next to it) are one
   finding. If points disagree, `point_mismatches`.
3. Per family, check the category with the tools in the table, starting
   from the family's `examples` (each has `evt` and `s_range_mm`):
   `ray --evt EVT --smin .. --smax ..` and `debug_ray` with a range of at
   most 20 mm around the interval. For a navigation artifact, compare a
   failing entry with one where B enters the same kind of part correctly
   before concluding.
4. Measure each finding:
   - placement: `boundary_shifts --part NAME@ID` (or `NAME@*` for all
     copies). If it has fewer than about 10 shifted boundaries, get the
     part's placements with `part_info`, `rescan` aimed at the part
     (target = a placement's `bbox_centre_mm` and `bbox_half_mm`, a few
     thousand rays), then `boundary_shifts` on the rescan case
     (`out/NAME/investigation/RESCAN`). D is measured only along
     `determined_along`; say so when the faces seen do not fix all of D. A
     large translation residual means a size change or rotation (see
     `faces` and `rigid`), or copies moved differently (`per_placement`).
   - material: affected parts and their map material.
   - navigation artifact: the share of the mismatch and the solids involved
     (shape, thickness from `debug_ray` local bounding boxes).
   - share of mismatch: the family `share` values (sum the families of one
     finding).
5. `triage` only when families are large and their cause is still open:
   its rules and modes are correlations; confirm them with a replay.
6. `submit_diagnosis` with:
   ```json
   {"schema": "simcadval-diagnosis/1", "case": "out/NAME",
    "verdict": "one sentence",
    "findings": [{"category": "placement", "parts": ["component0@4"],
                  "materials": {"a": "Carbon", "b": "Air"},
                  "share_of_mismatch": 0.4,
                  "magnitude": {"summary": "D = (0, 0, +0.010) mm",
                                "displacement_mm": [0, 0, 0.010]},
                  "confidence": "high",
                  "explanation": "why this category, and what rules out the others",
                  "evidence": [{"call": "c004", "shows": "boundary_pattern one_end_each"}]}],
    "unexplained": {"share_of_mismatch": 0.0, "notes": ""},
    "next_steps": []}
   ```
   Parts are `name@id` (from the tools' output) or `name@*` for every
   part of that name. If the result has `ok: false`, fix every listed error
   and submit again. Then report the verdict and the path of
   `diagnosis.md`.

## Rules

- Every claim cites a call id whose output shows it; quote the values in
  `shows`.
- Do not report a cause you have not shown: use `unclear` and say what
  would settle it.
- At most 3 rescans per investigation, each 20000 rays or fewer.
- Keep G4OCCT internals out of the diagnosis beyond the observed answers of
  its solids; for a navigation artifact, point to eic/G4OCCT#405 when it
  matches, or say that it looks new.
