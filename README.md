# simcadval

Validate CAD (STEP) detector geometry against the simulation geometry it is
meant to replace. The two are compared by **the questions Geant4 asks of
them**, not by converting one into the other:

1. *"Where is the next boundary along this ray?"* (geantino tracking steps)
2. *"What material am I in?"* (`G4Navigator::LocateGlobalPointAndSetup`)

A DD4hep detector is loaded as a Geant4 world twice: once through DDG4's
`Geant4Converter` (what `ddsim` runs) and once from a STEP file through
[G4OCCT](https://github.com/eic/G4OCCT) (`G4OCCTAssemblyVolume`). The same
deterministic set of geantino rays and points goes through both worlds. If
the answers agree on a dense sample, the simulation cannot tell the two
worlds apart.

```
compact.xml ──dd4hep2step──▶ X.step + X.step.materials.xml + export.json
     │                                   │
     ▼                                   ▼
simcadval_scan --world dd4hep    simcadval_scan --world step    (same rays/points)
     │  dd4hep.{rays,steps,points}.bin   │  step.{rays,steps,points}.bin
     └───────────────┬───────────────────┘
                     ▼
               compare.py ──▶ metrics.json, summary.txt, plots, mismatch CSVs
```

## Requirements

DD4hep, Geant4, ROOT, OpenCASCADE, G4OCCT and Python 3 with numpy and
matplotlib. All of them ship in the EIC `eic_xl` container. `eic.sh` runs a
command inside it (needs singularity/apptainer and `/cvmfs`). The image is
pinned to `eic_xl:26.11.0-stable`; override it with `EIC_IMG=...`.

## Build

```bash
./eic.sh 'cmake -S . -B build && cmake --build build -j'
```

This builds against the G4OCCT that ships in the container. To use another
G4OCCT build, add `-DCMAKE_PREFIX_PATH=<its install prefix>`.

## Run

```bash
./eic.sh 'scripts/run_case.sh minisid examples/MiniSiD.xml \
          --rays 20000 --ray-mode ip --points 1000000 --box 1450,1450,1800'
```

`run_case.sh NAME COMPACT [scan options]` exports the STEP file, scans both
worlds and compares them. Results go to `out/NAME/`: `summary.txt`, plus
`cmp/metrics.json`, plots, and `ray_mismatches.csv` / `point_mismatches.csv`.

Scan options:

| Option | Meaning |
|---|---|
| `--rays N` | number of geantino rays |
| `--ray-mode ip` | isotropic from the origin |
| `--ray-mode box` | origin uniform in `--box`, isotropic direction |
| `--ray-mode aim` | origin in `--box`, aimed at a uniform point of `--target cx,cy,cz,hx,hy,hz` |
| `--points M` | points uniform in `--box hx,hy,hz` (mm) |
| `--seed S` | ray *i* and point *j* depend only on (seed, index), never on the world |

Environment overrides for `run_case.sh`: `STEP=<file>` (use an existing STEP
file instead of exporting one), `MATMAP=<file>` (material map for the STEP
world), `EXPORT_ARGS=...` (extra `dd4hep2step` options), `SCANNER=<exe>`,
`BUILD=<dir>`, `OUT=<dir>`.

To check the metrics themselves, `EXPORT_ARGS="--perturb <volume>:dx,dy,dz"`
exports deliberately misplaced CAD.

## Export (`dd4hep2step`)

* A volume becomes a STEP part if its material differs from the world
  material. World-material regions stay empty and resolve to the world on
  import.
* Each part is its own shape **minus all daughters** (assemblies looked
  through), so parts never overlap and G4OCCT's flat import is a valid
  placement.
* One prototype per `TGeoVolume`, placed by XCAF references (STEP instancing
  → shared logical volumes).
* Part names are `<volume>@<id>`. The XDE material attribute is the DD4hep
  material name. A G4OCCT `materials.xml` is written as well, because XDE
  materials do not survive the OCCT 8.0 write→read round trip.
* Shapes: Box, Tube/TubeSeg/Cone/ConeSeg/Pcon, Pgon, Trd1/Trd2/Arb8/Trap/Gtra/Para,
  Xtru, Sphere, Torus, Eltu, Boolean composites, Scaled.
* Every part is read back from the written STEP file. Each one must come back
  as a single solid whose volume matches `TGeoShape::Capacity()`. If any does
  not, the exporter exits with code 4.

## Metrics (`compare.py`)

Ray metrics overlay the two piecewise-constant material functions along each
ray exactly (the union of both sets of step boundaries), so no along-ray
sampling step enters. "Non-world length" is the ray length where at least one
world is not in the world material. It is the denominator for the fractions,
so air outside the detector does not dilute them.

| Metric | Definition |
|---|---|
| material agreement | 1 − (length with mat_A ≠ mat_B) / non-world length |
| per-material IoU | L(A=m ∧ B=m) / L(A=m ∨ B=m) |
| confusion matrix | ray length with (mat_A, mat_B) |
| X/X₀, λ/λ_I per ray | Σ len/X₀ in each world, and their relative difference |
| boundary P/R/F1 @ τ | material-change boundaries matched to the nearest same-transition boundary on the same ray within τ (1 nm … 1 mm) |
| max boundary shift | largest nearest-boundary distance along rays. It is a real shift only when the transitions match; for a missing boundary it is the gap to the next one |
| sequence identity | fraction of rays with an identical collapsed material sequence |
| volume-id agreement | material agreement computed on DD4hep volume ids instead of materials |
| point agreement | fraction of points with equal material (Wilson 95 % CI), and its confusion matrix |
| cost | CPU µs per ray and per point in each world |

## Layout

| Path | What |
|---|---|
| `src/dd4hep2step.cc`, `tgeo_reader.cc`, `occt_writer.cc`, `geo_desc.hh` | DD4hep → STEP exporter. ROOT and OCCT live in separate translation units (`geo_desc.hh` is the std-only interface between them) |
| `src/simcadval_scan.cc`, `step_world.cc` | Geant4 app: one world per run, geantino scan, point sampling, `--debug-ray` replay |
| `src/geo_ids.hh` | volume numbering shared by the exporter and the scanner |
| `scripts/run_case.sh` | export → scan ×2 → compare, for one case |
| `scripts/compare.py` | metrics and plots for one case |
| `examples/MiniSiD.xml` | trimmed DD4hep SiD: beampipe, Si vertex barrel, 12-stave ECal barrel |
| `eic.sh` | runs a command in the `eic_xl` container |
