#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""check_constant_case.py -- a named constant is spelled by its kind, not by who uses it.

  * The physical, mathematical and unit-conversion constants and the kinds, in meds_constants and
    meds_kinds, are lower case: they read as symbols in a formula (grav, rho_h2o, day_sec, wp).
  * Every other named constant is UPPER case: a selector code, an index or a flag (SOIL_RETENTION_VG,
    INTEG_ARK), and a number chosen in the code -- a tolerance, a floor, a table size or a fitted
    coefficient (RWC_FLOOR, HYDRO_TABLE_NTAB) -- public or private, module-level or local.

A constant keeps its name when another module starts to use it. Fortran ignores case, so the
spelling is for the reader only, and this check keeps it one way.

Usage: check_constant_case.py REPO_ROOT   (exit 1 and a list when a constant breaks the rule)
"""
import re
import sys
from pathlib import Path

LOWER_FILES = {"src/share/meds_constants.f90", "src/share/meds_kinds.f90"}
DECL = re.compile(r"^\s*(?:real|integer|logical|character|complex|type)\b(.*?)::(.*)$", re.I)


def strip_comment(line):
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            return line[:i]
    return line


def names(entities):
    """The names declared in `a = 1, b(3) = [...], c = 'x, y'`: split on top-level commas."""
    out, depth, cur, quote = [], 0, "", None
    for ch in entities + ",":
        if quote:
            cur += ch
            quote = None if ch == quote else quote
            continue
        if ch in "'\"":
            quote = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            m = re.match(r"\s*(\w+)", cur)
            if m:
                out.append(m.group(1))
            cur = ""
        else:
            cur += ch
    return out


def main(root):
    bad = []
    files = sorted(list(root.glob("src/**/*.f90")) + list(root.glob("src/**/*.F90")) + list(root.glob("test/**/*.f90")))
    for path in files:
        rel = path.relative_to(root).as_posix()
        lines = path.read_text().split("\n")
        i = 0
        while i < len(lines):
            code, j = strip_comment(lines[i]), i
            while code.rstrip().endswith("&") and j + 1 < len(lines):
                j += 1
                code = code.rstrip()[:-1] + " " + strip_comment(lines[j]).lstrip().lstrip("&")
            m = DECL.match(code)
            if m and re.search(r"\bparameter\b", m.group(1), re.I):
                for name in names(m.group(2)):
                    want = name.lower() if rel in LOWER_FILES else name.upper()
                    if name != want:
                        bad.append(f"{rel}:{i + 1}: {name} -> {want}")
            i = j + 1
    if bad:
        print(f"{len(bad)} named constant(s) break the case rule (lower case only in meds_constants and "
              "meds_kinds, UPPER everywhere else):")
        print("\n".join("  " + b for b in bad))
        return 1
    print(f"constants named by kind: {len(files)} files checked")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]).resolve()))
