# SPDX-License-Identifier: Apache-2.0
"""TOML in and out for trial configs.

A trial's configs are made by PARSING the base TOML and re-serializing it, never by patching text
(scripts/numerics_sweep.py does the same, and says why: an appended block that repeats a table
re-parents the keys after it). MEDS reads a small TOML subset -- `key = value` lines under
`[section]` or `[section.sub]` headers, with scalars, strings and flat arrays -- and this writes
exactly that.
"""
from __future__ import annotations

import copy
import datetime as _dt

try:
    import tomllib  # py3.11+
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib


def load(path) -> dict:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    if isinstance(v, (_dt.datetime, _dt.date, _dt.time)):
        return v.isoformat()
    raise TypeError(f"cannot serialize {type(v).__name__} to TOML")


def dumps(d: dict) -> str:
    """Nested dicts to TOML, each table's scalars before its sub-tables."""
    def walk(table: dict, prefix: str, out: list):
        scalars = {k: v for k, v in table.items() if not isinstance(v, dict)}
        tables = {k: v for k, v in table.items() if isinstance(v, dict)}
        if prefix and (scalars or not tables):
            out.append(f"[{prefix}]")
        for k, v in scalars.items():
            out.append(f"{k} = {_fmt(v)}")
        if scalars:
            out.append("")
        for k, v in tables.items():
            walk(v, f"{prefix}.{k}" if prefix else k, out)
    lines: list[str] = []
    walk(d, "", lines)
    return "\n".join(lines) + "\n"


def write(path, d: dict) -> None:
    with open(path, "w") as fh:
        fh.write(dumps(d))


def deep_get(d: dict, key: str, default=None):
    node = d
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


def deep_set(d: dict, key: str, value) -> None:
    node = d
    parts = key.split(".")
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    node[parts[-1]] = value


def clone(d: dict) -> dict:
    return copy.deepcopy(d)


def digest(*objs) -> str:
    """A content hash that ignores key order: the same settings give the same digest however the
    tables were assembled."""
    import hashlib
    import json
    text = json.dumps(objs, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha1(text.encode()).hexdigest()
