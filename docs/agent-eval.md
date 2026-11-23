# Agent eval: runs so far, and why the skill and tools changed

The investigation procedure (`.claude/skills/investigate-mismatch/SKILL.md`)
and its tools (`scripts/svtool.py`) are written for any model. This page
records the blind runs on the one graded case so far, and which change each
run made necessary. Claude Code found everything (6/6); the best local model
so far scored 5/6 once, with the current skill.

## The graded case

`shift10um` from `scripts/selftest.sh`: in `examples/MiniSiD.xml`, the inner
vertex modules `component0@4` and `component1@5` are moved by
D = (0, 0, +0.010) mm in the STEP world. On top of that sits the known G4OCCT
navigation artifact (thin ECal slices read as Air, eic/G4OCCT#405), which is
99.9998 % of the mismatch length. A full answer names both, with the
moved parts and D.

The agent gets a stripped copy, because the selftest case directory gives
the answer away (`export.json` `perturb` entry, `export.log`,
`cmp_vs_base.txt`, the STEP file name). The copy has:

- a neutral directory name, and the STEP file and its `.materials.xml`
  renamed to match (`step.meta.json` rewritten to point at them);
- `export.json` without the `perturb` entry; no `export.log`, logs,
  `summary.txt`, `cmp_null*` or `cmp_vs_base*`;
- the scans hard-linked (they are large); every file a run writes is new,
  so the source case is not touched. Give each run its own copy.

Nothing stops an agent working in the repo from reading `scripts/grade.py`
or `scripts/selftest.sh`, which hold the answer; check the call log
(`investigation/calls.jsonl`) or the transcript for that.

## Grading

```
scripts/grade.py out/eval/NAME --answer shift10um
```

Six checks on `cmp/diagnosis.json`: navigation artifact found (share >= 0.9);
its parts are ECal slices; placement found; placement parts are only the
moved modules; both modules named; D within 0.002 mm in z and 0.003 mm in x
and y. Minus one per false finding (another category with share > 0.01).

## Local models

```
scripts/ollama_agent.py --model MODEL --case out/eval/NAME --answer shift10um
```

The same procedure and tools as the skill: the system prompt is the
SKILL.md body, the tools come from `svtool.py spec --format openai`, and
each call runs `svtool.py call`. It needs an Ollama server (default
`http://127.0.0.1:11434`), allows at most 3 rescans, and writes
`investigation/agent_<model>.jsonl`. `--skill FILE` runs another version of
the procedure.

## Runs

All on 2026-11-16, Ollama 0.33.3, `--num-ctx 32768`. "Skill" is the version
of SKILL.md: the svtool rewrite, or that plus the rounding rule below.

| Harness, model | Skill | Score /6 | Calls | Minutes | What happened |
|---|---|---|---|---|---|
| Claude Code, claude-opus-5-5 | rewrite | 6 | 24 (2 rescans) | ~19 | everything, D = (0, 0, +0.010) mm |
| Ollama, Qwen3.6-35B-A3B-MXFP4 | rewrite | 2 | 9 | 2.2 | artifact right; called the 0.01 mm families "floating point" |
| Ollama, Qwen3.6-35B-A3B-MXFP4 | + rounding | 5 | 13 | 1.7 | artifact and D = (0, 0, 0.01) mm; named `component0@4` only |
| Ollama, Qwen3.6-35B-A3B-MXFP4 | + rounding | 2 | 7 | 1.3 | stopped after the artifact |
| Ollama, Nemotron-Orchestrator-8B Q4_K_M | + rounding | 2 | 8 | 5.2 | artifact right; called the shift an extra part (share below the false-finding cut) |
| Ollama, qwen3.5:4b | rewrite | 1 | 14 | 12.4 | stopped at the first family |
| Ollama, gpt-oss:20b (think medium) | + rounding | -1 | 5 | 0.9 | a false missing part |
| Ollama, gpt-oss:20b (think medium) | rewrite | -2 | 23 | 3.5 | 13 `debug_ray` calls with `world` A/B came back empty; artifact called "material" |
| Ollama, gpt-oss:20b (think medium), 2 runs | rewrite | none | 5 each | 6-7 | Ollama answered HTTP 500 on its tool calls; no diagnosis |
| Ollama, Nemotron-Orchestrator-8B Q4_K_M | rewrite | none | 4 | 53 | Ollama put the model mostly on CPU; no diagnosis (not a fair test) |

Before the rewrite, the first skill (shell commands, hand-written
`diagnosis.md`) found both faults in Claude Code with 44 tool calls, about
127k tokens and 11 minutes. That run had the unstripped case and is not
graded.

## Why each change was needed

| Change | Needed because |
|---|---|
| Skill drives svtool and submits a checked diagnosis | The first skill had the agent run shell commands and dig small families out of `features.csv.gz`. A local model cannot do that through function calls, and a hand-written `diagnosis.md` cannot be graded. With the rewrite, Claude Code scored 6/6 (not faster: 24 calls and ~19 min against 44 and 11), and every harness now produces `diagnosis.json` that `grade.py` scores. |
| `svtool` rejects values outside an argument's enum, and maps `A`/`B` to `dd4hep`/`step` | gpt-oss passed `world: "B"` to `debug_ray`. The scanner got `--world B`, found no points and reported none, without an error, so the model lost 13 calls. Now `A`/`B` work and anything else is an error the model can read. |
| Skill: offsets of 1e-6 mm or more are real; `one_end_each` is never "precision" | Qwen3.6 found the 0.01 mm families but dismissed them as floating point, because their share of the mismatch is tiny (about 1.5e-6). The tools already drop rounding (below 1e-7 mm), so a family they report is geometry. Without the rule the one Qwen3.6 run missed the shift; with it, one of two runs measured it. |
| `ollama_agent.py` and `grade.py` in the repo | They produced and scored the runs above; without them the numbers cannot be reproduced. |

## Open

- One case only. The two Qwen3.6 runs with the rounding rule scored 5 and 2,
  so a score needs several runs per model.
- Claude Code has not been run on the skill with the rounding rule.
- The other planned eval cases (rotation, scale, missing part, material
  swap) need new `dd4hep2step` options first.
