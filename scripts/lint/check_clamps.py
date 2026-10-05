#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""check_clamps.py -- no bare literal clamps in the model source (#346).

A max/min against a literal hides a threshold in the arithmetic. The rule:
  * a guard against dividing by zero or leaving a function's domain uses tiny_num;
  * a physical threshold is a named setting or parameter with its source;
  * 0 and +-1 are physical bounds (non-negativity, fractions, the domain of acos) and pass.
A literal that is right where it is (a solver's step size, a function's domain) carries the
reason on the same line: `! clamp-ok: <reason>`.

Usage: check_clamps.py SRC_DIR   (exit 1 and a list when a clamp is bare)
"""
import re
import sys
from pathlib import Path

NUM = r"[-+]?(?:\d+\.\d*(?:[eEdD][-+]?\d+)?|\.\d+(?:[eEdD][-+]?\d+)?|\d+[eEdD][-+]?\d+)(?:_\w+)?"   # reals only
CALL = re.compile(r"\b(max|min)\s*\(", re.I)
ALLOWED = re.compile(r"^[-+]?(?:0|1)(?:\.0*)?(?:[eEdD][-+]?0+)?(?:_\w+)?$")


def args_of(text, start):
    """The top-level arguments of the call whose '(' is at text[start]; None if it never closes."""
    depth, cur, out = 0, "", []
    for ch in text[start:]:
        if ch == "(":
            depth += 1
            if depth == 1:
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                out.append(cur.strip())
                return out
        if depth == 1 and ch == ",":
            out.append(cur.strip())
            cur = ""
        elif depth >= 1:
            cur += ch
    return None


def logical_lines(path):
    """Fortran statements joined across '&' continuations, with their first line number."""
    buf, first = "", None
    for no, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        code, _, comment = raw.partition("!")
        if first is None:
            first = no
        tag = "clamp-ok:" in comment
        buf += " " + code.strip().lstrip("&") + (" !clamp-ok" if tag else "")
        if code.rstrip().endswith("&"):
            buf = buf.rstrip().rstrip("&")
            continue
        yield first, buf
        buf, first = "", None


def bare_clamps(src: Path):
    found = []
    for f in sorted(src.rglob("*.f90")):
        for no, stmt in logical_lines(f):
            if "!clamp-ok" in stmt:
                continue
            for m in CALL.finditer(stmt):
                args = args_of(stmt, m.end() - 1)
                if not args:
                    continue
                for a in args:
                    if re.fullmatch(NUM, a) and not ALLOWED.match(a):
                        found.append(f"{f}:{no}: {m.group(1)}(..., {a}): {stmt.strip()[:110]}")
    return found


if __name__ == "__main__":
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "src")
    bad = bare_clamps(src)
    for line in bad:
        print(line)
    if bad:
        print(f"\n{len(bad)} bare literal clamp(s): use tiny_num for a guard, a named setting for a "
              "threshold, or say why on the line with '! clamp-ok: <reason>' (scripts/lint/check_clamps.py)")
        sys.exit(1)
    print("check_clamps: no bare literal clamps")
