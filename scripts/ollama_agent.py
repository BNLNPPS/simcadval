#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run the investigate-mismatch procedure with a local model served by Ollama (stdlib only).

    ollama_agent.py --model MODEL --case out/eval/NAME [--answer shift10um] [--num-ctx 32768]

The same procedure and tools as the Claude Code skill, so runs can be compared: system prompt =
the SKILL.md body, tools = `svtool.py spec --format openai`, every tool call runs
`scripts/svtool.py call NAME JSON`. Writes CASE/investigation/agent_<model>.jsonl and prints a JSON
summary; with --answer it also grades the diagnosis (scripts/grade.py). See docs/agent-eval.md.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SVTOOL = os.path.join(HERE, "svtool.py")
sys.path.insert(0, HERE)
from grade import ANSWERS, grade  # noqa: E402
PREFACE = """You are an investigation agent for the simcadval project. You act only by calling the provided \
tools (function calls); where the procedure below shows shell commands like `scripts/svtool.py TOOL --arg value`, \
call the tool TOOL with those arguments instead. Work on case {case}: pass "case": "{case}" to every tool. \
Finish by calling submit_diagnosis; if it returns ok:false, fix every error and call it again. When it returns \
ok:true, reply with the verdict in one sentence and stop.

"""


def chat(url, body, timeout=1800):
    req = urllib.request.Request(url + "/api/chat", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def run_tool(name, args):
    try:
        p = subprocess.run([SVTOOL, "call", name, json.dumps(args)], capture_output=True, text=True,
                           timeout=3600, cwd=REPO)
        out = p.stdout.strip().splitlines()
        return out[-1] if out else json.dumps({"error": f"no output (rc {p.returncode}): {p.stderr[-300:]}"})
    except subprocess.TimeoutExpired:
        return json.dumps({"error": "tool timed out"})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--url", default="http://127.0.0.1:11434")
    ap.add_argument("--num-ctx", type=int, default=32768)
    ap.add_argument("--max-turns", type=int, default=60)
    ap.add_argument("--max-tool-chars", type=int, default=7000)
    ap.add_argument("--temperature", type=float, default=0.2)
    ap.add_argument("--think", default="true", help="true/false or low/medium/high (gpt-oss)")
    ap.add_argument("--wall-min", type=float, default=90)
    ap.add_argument("--skill", default=os.path.join(REPO, ".claude/skills/investigate-mismatch/SKILL.md"))
    ap.add_argument("--answer", choices=sorted(ANSWERS),
                    help="grade the diagnosis against this known answer (scripts/grade.py)")
    a = ap.parse_args()
    skill = open(a.skill).read()
    skill = skill.split("---", 2)[2] if skill.startswith("---") else skill
    tools = json.loads(subprocess.check_output([SVTOOL, "spec", "--format", "openai"], cwd=REPO))
    msgs = [{"role": "system", "content": PREFACE.format(case=a.case) + skill},
            {"role": "user", "content": f"Investigate the mismatch in {a.case}."}]
    think = {"true": True, "false": False}.get(a.think, a.think)
    budget = int(a.num_ctx * 2.6)  # chars of history kept in full, roughly 0.8 of num_ctx in tokens
    inv = os.path.join(REPO, a.case, "investigation")
    os.makedirs(inv, exist_ok=True)
    log = open(os.path.join(inv, f"agent_{re.sub(r'[^A-Za-z0-9.]+', '_', a.model)}.jsonl"), "w")
    t0, ok, rescans, nudges, ncalls, errors = time.time(), False, 0, 0, 0, 0
    tok_in = tok_out = 0
    final = ""
    for turn in range(a.max_turns):
        if time.time() - t0 > 60 * a.wall_min:
            break
        # keep the context bounded: shorten the oldest tool results first
        total = sum(len(m.get("content") or "") for m in msgs)
        for m in msgs:
            if total <= budget:
                break
            if m["role"] == "tool" and len(m["content"]) > 400:
                total -= len(m["content"]) - 400
                m["content"] = m["content"][:400] + " ...[older output shortened; call the tool again if needed]"
        body = {"model": a.model, "messages": msgs, "tools": tools, "stream": False, "think": think,
                "options": {"num_ctx": a.num_ctx, "temperature": a.temperature}}
        r = None
        for attempt in range(3):  # a server error (e.g. an unparsable tool call) is retried, warmer each time
            body["options"]["temperature"] = a.temperature + 0.4 * attempt
            try:
                r = chat(a.url, body)
                break
            except Exception as e:  # noqa: BLE001
                log.write(json.dumps({"turn": turn, "attempt": attempt, "error": str(e)}) + "\n")
                log.flush()
                time.sleep(5)
        if r is None:
            break
        m = r.get("message", {})
        tok_in += r.get("prompt_eval_count", 0)
        tok_out += r.get("eval_count", 0)
        calls = m.get("tool_calls") or []
        msgs.append({"role": "assistant", "content": m.get("content", ""), **({"tool_calls": calls} if calls else {})})
        rec = {"turn": turn, "t": round(time.time() - t0, 1), "prompt_tokens": r.get("prompt_eval_count"),
               "out_tokens": r.get("eval_count"), "content": (m.get("content") or "")[:600],
               "thinking": (m.get("thinking") or "")[:600], "calls": []}
        if not calls:
            log.write(json.dumps(rec) + "\n")
            log.flush()
            if ok:
                final = m.get("content", "")
                break
            nudges += 1
            if nudges > 4:
                break
            msgs.append({"role": "user", "content": "Continue the procedure by calling tools. No diagnosis has been "
                                                    "accepted yet: finish with submit_diagnosis."})
            continue
        for tc in calls:
            fn = tc.get("function", {})
            name, args = fn.get("name"), fn.get("arguments", {})
            if isinstance(args, str):
                try:
                    args = json.loads(args)
                except ValueError:
                    args = {}
            if name == "rescan" and rescans >= 3:
                res = json.dumps({"error": "rescan budget (3) used up"})
            else:
                rescans += name == "rescan"
                res = run_tool(name, args)
            ncalls += 1
            try:
                rj = json.loads(res)
                errors += "error" in rj
                if name == "submit_diagnosis" and rj.get("ok"):
                    ok = True
            except ValueError:
                errors += 1
            if len(res) > a.max_tool_chars:
                res = res[:a.max_tool_chars] + " ...[truncated]"
            msgs.append({"role": "tool", "content": res, "tool_name": name})
            rec["calls"].append({"name": name, "args": args, "result": res[:500]})
        log.write(json.dumps(rec) + "\n")
        log.flush()
    summary = {"model": a.model, "case": a.case, "accepted": ok, "turns": turn + 1, "tool_calls": ncalls,
               "tool_errors": errors, "rescans": rescans, "nudges": nudges, "prompt_tokens": tok_in,
               "output_tokens": tok_out, "minutes": round((time.time() - t0) / 60, 1), "final": final[:300]}
    if a.answer:
        summary["grade"] = grade(os.path.join(REPO, a.case), a.answer)
    log.write(json.dumps({"summary": summary}) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
