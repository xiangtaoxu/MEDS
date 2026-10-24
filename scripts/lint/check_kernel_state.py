#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""check_kernel_state.py -- the kernels never see polygon_t.

A kernel works on one patch at a time, through the store types that describe what one patch holds
(canopy air, snow, soil water and heat, soil carbon, and the two accumulators). It never sees the
polygon, the cohort block or the patch block: that keeps a kernel testable alone and eligible for an
OpenMP `target` device. A routine that needs more is driver code.

The library boundary cannot hold this rule, because every module file goes to one shared
directory: once the state layer is built, any file can `use` any of its modules. So this script
checks it on every ctest run. Code in the kernel folders (src/fast_dynamics/processes,
src/slow_dynamics/processes):
  * takes from the state layer only the per-patch store types and their blends, through
    `use meds_patch_state_vars, only : ...`;
  * never names polygon_t, cohort_block or patch_block;
  * never uses meds_polygon_state_vars, meds_cohort_state_vars or meds_diag_vars.

Usage: check_kernel_state.py SRC_DIR   (exit 1 and a list when a kernel reaches past its patch)
"""
import re
import sys
from pathlib import Path

KERNEL_DIRS = ("fast_dynamics/processes", "slow_dynamics/processes")
FORBIDDEN_TYPES = ("polygon_t", "cohort_block", "patch_block")
FORBIDDEN_MODULES = ("meds_polygon_state_vars", "meds_cohort_state_vars", "meds_diag_vars")
PATCH_MODULE = "meds_patch_state_vars"
ALLOWED = {"cas_state_t", "snow_state_t", "soil_water_state_t", "soil_energy_state_t",
           "soil_carbon_state_t",
           "decomp_xi_t", "litter_input_t", "blend_cas", "blend_snow", "blend_soil_water", "blend_soil_energy",
           "blend_soil_carbon", "blend_decomp_xi", "blend_litter_input", "cas_set_depth",
           "cas_capacities", "soil_water_store", "soil_energy_store", "pack_pool_vector", "unpack_pool_vector"}
USE = re.compile(r"^\s*use\s*(?:,\s*\w+\s*::\s*|::\s*|\s)\s*(\w+)\s*(,\s*only\s*:(.*))?$", re.I)


def strip_comment(line):
    """The line without its `!` comment (a `!` inside a quoted string is kept)."""
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


def statements(path):
    """(first line number, statement) for each statement, continuation lines joined."""
    buf, start = "", 0
    for no, raw in enumerate(path.read_text(errors="replace").split("\n"), 1):
        code = strip_comment(raw).rstrip()
        if not buf:
            start = no
        cont = code.endswith("&")
        buf += " " + (code[:-1] if cont else code).strip().lstrip("&")
        if not cont:
            if buf.strip():
                yield start, buf.strip()
            buf = ""


def violations(src):
    found = []
    for d in KERNEL_DIRS:
        for f in sorted((src / d).glob("*.f90")):
            for no, stmt in statements(f):
                for t in FORBIDDEN_TYPES:
                    if re.search(rf"\b{t}\b", stmt, re.I):
                        found.append(f"{f}:{no}: names {t}: {stmt[:100]}")
                m = USE.match(stmt)
                if not m:
                    continue
                mod = m.group(1).lower()
                if mod in FORBIDDEN_MODULES:
                    found.append(f"{f}:{no}: uses {mod}")
                elif mod == PATCH_MODULE:
                    if not m.group(2):
                        found.append(f"{f}:{no}: uses {mod} without an only-list")
                        continue
                    names = [n.strip().split("=>")[-1].strip().lower() for n in m.group(3).split(",") if n.strip()]
                    for n in names:
                        if n not in ALLOWED:
                            found.append(f"{f}:{no}: takes {n} from {mod} (only the per-patch store types are allowed)")
    return found


if __name__ == "__main__":
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "src")
    bad = violations(src)
    for line in bad:
        print(line)
    if bad:
        print(f"\n{len(bad)} kernel reference(s) past the patch: a kernel takes only the per-patch store "
              "types; a routine that needs more is driver code (scripts/lint/check_kernel_state.py)")
        sys.exit(1)
    print("check_kernel_state: the kernels take only the per-patch store types")
